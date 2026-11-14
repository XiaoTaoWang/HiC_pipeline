#!/usr/bin/env python3
"""Bounded-memory, fixed-bin Pairix -> Cooler loader.

Requires Python >=3.10, cooler, pypairix (and Cooler's numpy/pandas/h5py).
Input: sorted, BGZIP-compressed, Pairix-indexed point pairs, one copy per pair,
upper-triangular within each chromosome. Chromosome order may differ from the
chrom.sizes file. No deduplication/filtering/balancing is performed.

Example:
  python fast_pairix_to_cool.py --assembly hg38 --nproc 4 \
      hg38.chrom.sizes:10 in.pairs.gz out.10bp.cool

Cooler's public create_cooler API initializes the schema using a tiny template.
Large datasets are replaced/written in bounded blocks with h5py, avoiding the
full-genome pandas bin table required by the ordinary creation API. No package
files or private Cooler APIs are modified. Output count dtype is int64.
"""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor
import logging
import multiprocessing as mp
import os
from pathlib import Path
import tempfile
import time

# Avoid multiplying numerical-library threads across worker processes.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import cooler
import h5py
import numpy as np
import pandas as pd
import pypairix

VERSION = "1.0.0"
LOG = logging.getLogger("fast_pairix")
RUN_DTYPE = np.dtype([("key", "<u8"), ("count", "<i8")])
_CONFIG = None
_PAIRIX = None


def read_bins_spec(spec):
    try:
        filename, size = spec.rsplit(":", 1)
        size = int(size)
    except ValueError as exc:
        raise ValueError("BINS must be a chrom.sizes filename followed by :BIN_SIZE") from exc
    if size <= 0:
        raise ValueError("Bin size must be positive")
    names, lengths = [], []
    with open(filename) as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(f"Expected two columns in chrom.sizes: {line.strip()}")
            name, length = fields[0], int(fields[1])
            if name in names or length <= 0 or length > np.iinfo(np.int32).max:
                raise ValueError(f"Duplicate chromosome or unsupported length: {name} {length}")
            names.append(name)
            lengths.append(length)
    if not names:
        raise ValueError("Empty chrom.sizes file")
    offsets = [0]
    for length in lengths:
        offsets.append(offsets[-1] + (length + size - 1) // size)
    # Packed unsigned keys are lossless if the full square fits in uint64.
    if offsets[-1] ** 2 > 2**64:
        raise ValueError("Too many bins for this script's uint64 pixel keys (maximum 2**32)")
    return names, lengths, size, offsets


def inspect_index(path, names, lengths, block_char):
    f = pypairix.open(path)
    columns = (f.get_startpos1_col(), f.get_startpos2_col())
    if (f.get_endpos1_col(), f.get_endpos2_col()) != columns:
        raise ValueError("Only point pairs are supported, not paired intervals")
    blocks = set()
    for block in f.get_blocknames():
        parts = block.split(block_char)
        if len(parts) != 2:
            raise ValueError("Cannot parse Pairix block names; check --block-char")
        blocks.add(tuple(parts))
    for a, b in blocks:
        if a != b and (b, a) in blocks:
            raise ValueError(f"Both {a}|{b} and {b}|{a} exist; input must be triangular")
    known = set(names)
    ignored = {c for block in blocks for c in block} - known
    if ignored:
        LOG.warning("Ignoring contacts involving chromosomes outside chrom.sizes: %s", ", ".join(sorted(ignored)))
    header_sizes = {}
    for line in f.get_header():
        if line.startswith("#chromsize:"):
            fields = line.split()
            if len(fields) == 3:
                header_sizes[fields[1]] = int(fields[2])
    for name, length in zip(names, lengths):
        if name in header_sizes and header_sizes[name] != length:
            raise ValueError(f"Header/chrom.sizes length mismatch for {name}")
    routes = []
    for i, a in enumerate(names):
        row = []
        for j in range(i, len(names)):
            b = names[j]
            if (a, b) in blocks:
                row.append((j, False))
            elif (b, a) in blocks:
                row.append((j, True))
        routes.append(row)
    del f
    return columns, routes


def init_worker(config):
    global _CONFIG, _PAIRIX
    _CONFIG = config
    _PAIRIX = pypairix.open(config["pairs"])


def combine_records(records):
    """Sort and sum a bounded batch of (packed key, count) records."""
    order = np.argsort(records["key"], kind="stable")
    records = records[order]
    starts = np.r_[0, np.flatnonzero(records["key"][1:] != records["key"][:-1]) + 1]
    result = np.empty(len(starts), dtype=RUN_DTYPE)
    result["key"] = records["key"][starts]
    result["count"] = np.add.reduceat(records["count"], starts)
    return result


def merge_runs(left, right, destination, block_size):
    """Merge two sorted unique runs, with bounded reads and vectorized merging.

    No Python loop per pixel and no mmap whose resident pages can accumulate.
    The smaller block endpoint is a safe cutoff: subsequent input keys exceed it.
    """
    with open(left, "rb") as lf, open(right, "rb") as rf, open(destination, "wb") as out:
        a = np.fromfile(lf, dtype=RUN_DTYPE, count=block_size)
        b = np.fromfile(rf, dtype=RUN_DTYPE, count=block_size)
        while len(a) and len(b):
            pivot = min(a["key"][-1], b["key"][-1])
            ia = np.searchsorted(a["key"], pivot, side="right")
            ib = np.searchsorted(b["key"], pivot, side="right")
            combine_records(np.concatenate((a[:ia], b[:ib]))).tofile(out)
            a, b = a[ia:], b[ib:]
            if not len(a):
                a = np.fromfile(lf, dtype=RUN_DTYPE, count=block_size)
            if not len(b):
                b = np.fromfile(rf, dtype=RUN_DTYPE, count=block_size)
        for pending, stream in ((a, lf), (b, rf)):
            while len(pending):
                pending.tofile(out)
                pending = np.fromfile(stream, dtype=RUN_DTYPE, count=block_size)
    os.unlink(left)
    os.unlink(right)


def process_interval(task):
    """Read all partners for one disjoint range of output matrix rows."""
    task_id, cid, lo, hi = task
    c, f = _CONFIG, _PAIRIX
    names, lengths, offsets = c["names"], c["lengths"], c["offsets"]
    size, nbins = c["binsize"], offsets[-1]
    p1, p2 = c["columns"]
    decrement = 0 if c["zero_based"] else 1
    keys = np.empty(c["buffer_contacts"], dtype=np.uint64)
    used = total = queries = serial = 0
    levels = []

    def new_path():
        nonlocal serial
        serial += 1
        return str(Path(c["temp_dir"]) / f"task-{task_id}-run-{serial}.bin")

    def flush(n):
        path = new_path()
        unique, counts = np.unique(keys[:n], return_counts=True)
        records = np.empty(len(unique), dtype=RUN_DTYPE)
        records["key"], records["count"] = unique, counts
        records.tofile(path)
        # Binary compaction: bounded metadata, and logarithmic merge depth.
        del unique, counts, records
        level = 0
        while level < len(levels) and levels[level] is not None:
            dest = new_path()
            merge_runs(levels[level], path, dest, c["io_chunk"])
            levels[level] = None
            path = dest
            level += 1
        if level == len(levels):
            levels.append(path)
        else:
            levels[level] = path

    for partner, flipped in c["routes"][cid]:
        queries += 1
        if flipped:
            iterator = f.query2D(names[partner], 0, lengths[partner], names[cid], lo, hi)
            row_col, col_col = p2, p1
        else:
            iterator = f.query2D(names[cid], lo, hi, names[partner], 0, lengths[partner])
            row_col, col_col = p1, p2
        for line in iterator:
            row_pos = int(line[row_col]) - decrement
            col_pos = int(line[col_col]) - decrement
            if not (lo <= row_pos < hi and 0 <= col_pos < lengths[partner]):
                raise ValueError("Out-of-range pair or index/coordinate mismatch; check chrom.sizes and --zero-based")
            if cid == partner and row_pos > col_pos:
                raise ValueError("Cis pairs must be upper-triangular (pos1 <= pos2)")
            i = offsets[cid] + row_pos // size
            j = offsets[partner] + col_pos // size
            keys[used] = i * nbins + j
            used += 1
            total += 1
            if used == len(keys):
                flush(used)
                used = 0
    if used:
        flush(used)
    del keys
    runs = [path for path in levels if path is not None]
    if not runs:
        return None, 0, queries
    path = runs[0]
    for other in runs[1:]:
        dest = new_path()
        merge_runs(path, other, dest, c["io_chunk"])
        path = dest
    return path, total, queries


def tasks_for(config, window_bp):
    # Every boundary is a bin boundary, so a pixel never spans two tasks.
    size = config["binsize"]
    width = max(1, window_bp // size) * size
    task_id = 0
    for cid, length in enumerate(config["lengths"]):
        if not config["routes"][cid]:
            continue
        for lo in range(0, length, width):
            yield task_id, cid, lo, min(lo + width, length)
            task_id += 1


def ordered_results(config, tasks, nproc):
    if nproc == 1:
        init_worker(config)
        for task in tasks:
            yield task, process_interval(task)
        return
    # Spawn avoids inheriting an open HDF5 writer. Bound submitted tasks, not
    # just worker count: later finished results must not accumulate unboundedly.
    with ProcessPoolExecutor(max_workers=nproc, mp_context=mp.get_context("spawn"),
                             initializer=init_worker, initargs=(config,)) as pool:
        pending = deque()
        tasks = iter(tasks)
        for _ in range(2 * nproc):
            task = next(tasks, None)
            if task is None:
                break
            pending.append((task, pool.submit(process_interval, task)))
        while pending:
            task, future = pending.popleft()
            yield task, future.result()
            task = next(tasks, None)
            if task is not None:
                pending.append((task, pool.submit(process_interval, task)))


def initialize_output(path, config, assembly, compression):
    """Use public Cooler API for schema, then stream the large fixed bin table."""
    names, lengths = config["names"], config["lengths"]
    seed = pd.DataFrame({"chrom": pd.Categorical(names, categories=names, ordered=True),
                         "start": np.zeros(len(names), dtype=np.int64), "end": lengths})
    empty = pd.DataFrame({name: pd.Series(dtype="int64")
                          for name in ("bin1_id", "bin2_id", "count")})
    cooler.create_cooler(path, seed, empty, ordered=True, assembly=assembly,
                         dtypes={"count": np.int64},
                         metadata={"loader": "fast_pairix_to_cool", "loader-version": VERSION},
                         h5opts={"compression": "gzip", "compression_opts": compression, "shuffle": True})
    nbins, block = config["offsets"][-1], config["io_chunk"]
    options = {"compression": "gzip", "compression_opts": compression, "shuffle": True,
               "chunks": (min(block, nbins),)}
    with h5py.File(path, "r+") as f:
        chrom_dtype = f["bins/chrom"].dtype
        for field, dtype in (("chrom", chrom_dtype), ("start", "i4"), ("end", "i4")):
            del f[f"bins/{field}"]
            f["bins"].create_dataset(field, shape=(nbins,), dtype=dtype, **options)
        for cid, length in enumerate(lengths):
            first, stop = config["offsets"][cid:cid + 2]
            for a in range(first, stop, block):
                b = min(a + block, stop)
                starts = (np.arange(a - first, b - first, dtype=np.int64) * config["binsize"])
                f["bins/chrom"][a:b] = np.full(b - a, cid, dtype=np.int32)
                f["bins/start"][a:b] = starts
                f["bins/end"][a:b] = np.minimum(starts + config["binsize"], length)
        del f["indexes/bin1_offset"]
        f["indexes"].create_dataset("bin1_offset", shape=(nbins + 1,), dtype="i8", **options)
        f["indexes/chrom_offset"][:] = config["offsets"]
        for field in ("bin1_id", "bin2_id", "count"):
            del f[f"pixels/{field}"]
            f["pixels"].create_dataset(field, shape=(0,), maxshape=(None,), dtype="i8",
                                        chunks=(block,), compression="gzip",
                                        compression_opts=compression, shuffle=True)
        f.attrs["nbins"] = nbins
        f.attrs["bin-type"] = "fixed"
        f.attrs["bin-size"] = config["binsize"]


def build_row_index(f, nbins, block):
    """Generate CSR row offsets in bounded blocks, including every empty row."""
    ids = f["pixels/bin1_id"]
    out = f["indexes/bin1_offset"]
    next_row = 0
    for pixel_start in range(0, len(ids), block):
        values = ids[pixel_start:pixel_start + block]
        last_row = int(values[-1])
        # The first occurrence of each newly seen row is in this pixel block.
        # Rows spanning pixel blocks were already assigned their first offset.
        for a in range(next_row, last_row + 1, block):
            b = min(a + block, last_row + 1)
            out[a:b] = pixel_start + np.searchsorted(values, np.arange(a, b), side="left")
        next_row = max(next_row, last_row + 1)
    for a in range(next_row, nbins + 1, block):
        b = min(a + block, nbins + 1)
        out[a:b] = np.full(b - a, len(ids), dtype=np.int64)


def write_pixels(path, config, results):
    nbins, block = config["offsets"][-1], config["io_chunk"]
    nnz = total = query_total = written_total = 0
    previous_key = -1
    last_log = time.monotonic()
    with h5py.File(path, "r+") as f:
        datasets = [f[f"pixels/{name}"] for name in ("bin1_id", "bin2_id", "count")]
        for task, (run, ncontacts, queries) in results:
            total += ncontacts
            query_total += queries
            if total > np.iinfo(np.int64).max:
                raise OverflowError("Contact total exceeds int64")
            if run is not None:
                with open(run, "rb") as stream:
                    while True:
                        records = np.fromfile(stream, dtype=RUN_DTYPE, count=block)
                        if not len(records):
                            break
                        keys = records["key"]
                        if int(keys[0]) <= previous_key or np.any(keys[1:] <= keys[:-1]):
                            raise RuntimeError("Internal error: pixels are not unique and globally sorted")
                        previous_key = int(keys[-1])
                        written_total += int(records["count"].sum())
                        arrays = ((keys // np.uint64(nbins)).astype(np.int64),
                                  (keys % np.uint64(nbins)).astype(np.int64), records["count"])
                        for dataset, values in zip(datasets, arrays):
                            dataset.resize((nnz + len(records),))
                            dataset[nnz:nnz + len(records)] = values
                        nnz += len(records)
                os.unlink(run)
            if time.monotonic() - last_log >= 10 or task[0] == 0:
                LOG.info("Completed interval %d (%s:%d-%d): %s contacts, %s pixels so far",
                         task[0] + 1, config["names"][task[1]], task[2], task[3], f"{total:,}", f"{nnz:,}")
                last_log = time.monotonic()
        if written_total != total:
            raise RuntimeError("Internal error: output counts do not match selected input contacts")
        LOG.info("Writing %s row offsets in blocks", f"{nbins + 1:,}")
        build_row_index(f, nbins, block)
        f.attrs["nnz"], f.attrs["sum"] = nnz, total
    return nnz, total, query_total


def run(args):
    names, lengths, binsize, offsets = read_bins_spec(args.bins)
    output = Path(args.output).absolute()
    pairs = str(Path(args.pairs).resolve())
    if "::" in args.output:
        raise ValueError("Output must be a standalone .cool filename, not an HDF5 group URI")
    if output.resolve() == Path(pairs) or output.resolve() == Path(args.bins.rsplit(":", 1)[0]).resolve():
        raise ValueError("Output must not overwrite an input")
    if output.exists() and not args.force:
        raise FileExistsError(f"Output already exists: {output}; use --force to replace")
    columns, routes = inspect_index(pairs, names, lengths, args.block_char)
    if args.temp_dir:
        Path(args.temp_dir).mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    # Both output and scratch remain private until success. Existing output is
    # preserved if any processing step fails, even with --force.
    with tempfile.TemporaryDirectory(prefix="fast-pairix-", dir=args.temp_dir) as scratch:
        config = {"names": names, "lengths": lengths, "binsize": binsize, "offsets": offsets,
                  "columns": columns, "routes": routes, "pairs": pairs, "temp_dir": scratch,
                  "zero_based": args.zero_based, "buffer_contacts": args.buffer_contacts,
                  "io_chunk": args.io_chunk}
        fd, partial = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".partial", dir=output.parent)
        os.close(fd)
        try:
            LOG.info("Preparing %s bins; %d workers; %s contacts per worker buffer",
                     f"{offsets[-1]:,}", args.nproc, f"{args.buffer_contacts:,}")
            initialize_output(partial, config, args.assembly, args.compression)
            results = ordered_results(config, tasks_for(config, args.window_bp), args.nproc)
            try:
                nnz, total, queries = write_pixels(partial, config, results)
            finally:
                results.close()
            check = cooler.Cooler(partial)
            if (not cooler.fileops.is_cooler(partial) or check.info["nbins"] != offsets[-1]
                    or check.info["nnz"] != nnz or check.binsize != binsize):
                raise RuntimeError("Cooler compatibility check failed")
            if args.force:
                os.replace(partial, output)
            else:
                # Same directory/filesystem: atomic publish that cannot clobber
                # an output created by another process while this job ran.
                os.link(partial, output)
                os.unlink(partial)
            LOG.info("Wrote %s: %s contacts, %s pixels, %s indexed queries in %.1f minutes",
                     output, f"{total:,}", f"{nnz:,}", f"{queries:,}", (time.monotonic() - start) / 60)
        finally:
            if os.path.exists(partial):
                os.unlink(partial)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bins", help="chrom.sizes:BIN_SIZE, e.g. hg38.chrom.sizes:10")
    parser.add_argument("pairs", help="Sorted BGZIP pairs with existing .px2 index")
    parser.add_argument("output", help="Standalone output .cool")
    parser.add_argument("--assembly", default=None)
    parser.add_argument("--nproc", type=int, default=4)
    parser.add_argument("--window-bp", type=int, default=1_000_000,
                        help="Query interval width, rounded down to a bin multiple (default 1000000)")
    parser.add_argument("--buffer-contacts", type=int, default=2_000_000,
                        help="Maximum packed contacts per worker buffer (default 2000000)")
    parser.add_argument("--io-chunk", type=int, default=250_000,
                        help="Bounded merge/index/output batch size (default 250000)")
    parser.add_argument("--temp-dir", help="Scratch directory, preferably on local SSD; default system temp")
    parser.add_argument("--compression", type=int, choices=range(10), default=1,
                        help="Output gzip level 0-9 (default 1)")
    parser.add_argument("--zero-based", action="store_true",
                        help="Input coordinates AND Pairix index are zero-based (default: 4DN one-based)")
    parser.add_argument("--block-char", default="|", help="Pairix index block separator")
    parser.add_argument("--force", action="store_true", help="Replace an existing output only after success")
    args = parser.parse_args(argv)
    for name in ("nproc", "window_bp", "buffer_contacts", "io_chunk"):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        run(args)
    except (ValueError, OSError, RuntimeError, pypairix.PairixError) as exc:
        LOG.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
