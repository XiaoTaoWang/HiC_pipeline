#!/usr/bin/env python3
"""Stream fixed-bin coarsening and perform in-memory sparse genome-wide ICE.

Python >=3.10; uses only cooler and its numpy/pandas/scipy/h5py dependencies.
Example:
  python fast_zoomify.py input.10bp.cool -o output.mcool --resolutions B \
      --threads 8 --memory-gb 160 --temp-dir /path/to/scratch

Every generated level, INCLUDING THE BASE, is balanced independently. Supports
full, fixed-width, symmetric-upper count matrices. No cis-only/trans-only mode.
Numerical ICE follows Cooler 0.10.4; floating-point reduction order differs.
Failed levels default to NaN weights and a nonzero final exit status. --resume
preserves completed raw levels and skips successfully balanced levels.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import gc
import json
import logging
import math
import os
from pathlib import Path
import tempfile
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import cooler
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

VERSION = "1.0.0"
LOG = logging.getLogger("fast_zoomify")
RECORD = np.dtype([("key", "<u8"), ("count", "<i8")])
GIB = 1024**3


def source_info(uri):
    c = cooler.Cooler(uri)
    info = c.info
    if not c.binsize or c.storage_mode != "symmetric-upper":
        raise ValueError("Input must have fixed-width bins and symmetric-upper storage")
    names = list(c.chromnames)
    lengths = c.chromsizes.to_numpy(dtype=np.int64)
    if len(lengths) == 0 or np.any(lengths <= 0) or np.any(lengths > np.iinfo(np.int32).max):
        raise ValueError("Chromosome lengths must be positive and fit int32")
    offsets = np.r_[0, np.cumsum((lengths + c.binsize - 1) // c.binsize)]
    if int(offsets[-1]) != info["nbins"] or int(offsets[-1]) > 2**32:
        raise ValueError("Input must cover complete chromosomes with at most 2**32 bins")
    with c.open() as g:
        if not np.array_equal(g["indexes/chrom_offset"][:], offsets):
            raise ValueError("Chromosome offsets do not match a complete fixed-bin grid")
        if g["pixels/count"].dtype.kind not in "iu":
            raise ValueError("This script supports nonnegative integer contact counts")
    return c, names, lengths, offsets


def resolutions_for(spec, base, genome_length):
    base = int(base)
    stop = math.ceil(genome_length / 256)
    if spec.upper() in ("B", "N"):
        values = [base]
        multiplier = 1
        while True:
            candidates = [base * multiplier] if spec.upper() == "B" else [base * multiplier * k for k in (1, 2, 5)]
            values.extend(r for r in candidates if base <= r <= stop)
            multiplier *= 2 if spec.upper() == "B" else 10
            if base * multiplier > stop:
                break
    else:
        values = [base] + [int(x.strip()) for x in spec.split(",")]
    values = sorted(set(values))
    if any(r < base or r % base for r in values):
        raise ValueError("Every resolution must be a positive integer multiple of the base")
    return values


def aggregate(records):
    if not len(records):
        return records
    order = np.argsort(records["key"], kind="stable")
    records = records[order]
    starts = np.r_[0, np.flatnonzero(records["key"][1:] != records["key"][:-1]) + 1]
    result = np.empty(len(starts), dtype=RECORD)
    result["key"] = records["key"][starts]
    result["count"] = np.add.reduceat(records["count"], starts)
    return result


def merge_runs(left, right, destination, block):
    with open(left, "rb") as lf, open(right, "rb") as rf, open(destination, "wb") as out:
        a = np.fromfile(lf, dtype=RECORD, count=block)
        b = np.fromfile(rf, dtype=RECORD, count=block)
        while len(a) and len(b):
            pivot = min(a["key"][-1], b["key"][-1])
            ia = np.searchsorted(a["key"], pivot, side="right")
            ib = np.searchsorted(b["key"], pivot, side="right")
            aggregate(np.concatenate((a[:ia], b[:ib]))).tofile(out)
            a, b = a[ia:], b[ib:]
            if not len(a):
                a = np.fromfile(lf, dtype=RECORD, count=block)
            if not len(b):
                b = np.fromfile(rf, dtype=RECORD, count=block)
        for data, stream in ((a, lf), (b, rf)):
            while len(data):
                data.tofile(out)
                data = np.fromfile(stream, dtype=RECORD, count=block)
    os.unlink(left)
    os.unlink(right)


def initialize_level(path, group, names, lengths, resolution, assembly, args):
    # Tiny public-API template; do not create a full-genome pandas bin table.
    seed = pd.DataFrame({"chrom": pd.Categorical(names, categories=names, ordered=True),
                         "start": np.zeros(len(names), dtype=np.int64), "end": lengths})
    empty = pd.DataFrame({k: pd.Series(dtype="int64") for k in ("bin1_id", "bin2_id", "count")})
    cooler.create_cooler(f"{path}::{group}", seed, empty, ordered=True, mode="a",
                         assembly=assembly, dtypes={"count": np.int64})
    offsets = np.r_[0, np.cumsum((lengths + resolution - 1) // resolution)]
    nbins, block = int(offsets[-1]), args.io_chunk
    options = dict(compression="gzip", compression_opts=args.compression, shuffle=True)
    with h5py.File(path, "r+") as f:
        g = f[group]
        chrom_dtype = g["bins/chrom"].dtype
        for name, dtype in (("chrom", chrom_dtype), ("start", "i4"), ("end", "i4")):
            del g[f"bins/{name}"]
            g["bins"].create_dataset(name, shape=(nbins,), dtype=dtype, chunks=(min(block, nbins),), **options)
        for cid, length in enumerate(lengths):
            first, end = map(int, offsets[cid:cid + 2])
            for a in range(first, end, block):
                b = min(a + block, end)
                starts = np.arange(a - first, b - first, dtype=np.int64) * resolution
                g["bins/chrom"][a:b] = np.full(b - a, cid, dtype=np.int32)
                g["bins/start"][a:b] = starts
                g["bins/end"][a:b] = np.minimum(starts + resolution, length)
        for name in ("bin1_id", "bin2_id", "count"):
            del g[f"pixels/{name}"]
            g["pixels"].create_dataset(name, shape=(0,), maxshape=(None,), dtype="i8", chunks=(block,), **options)
        del g["indexes/bin1_offset"]
        g["indexes"].create_dataset("bin1_offset", shape=(nbins + 1,), dtype="i8", chunks=(min(block, nbins + 1),), **options)
        g["indexes/chrom_offset"][:] = offsets
        g.attrs.update({"nbins": nbins, "bin-size": resolution, "bin-type": "fixed"})
    return offsets


def make_row_index(g, block):
    ids, offsets = g["pixels/bin1_id"], g["indexes/bin1_offset"]
    next_row = 0
    for start in range(0, len(ids), block):
        values = ids[start:start + block]
        end = int(values[-1]) + 1
        for a in range(next_row, end, block):
            b = min(a + block, end)
            offsets[a:b] = start + np.searchsorted(values, np.arange(a, b), side="left")
        next_row = max(next_row, end)
    for a in range(next_row, len(offsets), block):
        b = min(a + block, len(offsets))
        offsets[a:b] = np.full(b - a, len(ids), dtype=np.int64)


def coarsen_level(path, source, target, names, lengths, args):
    temporary_group = "/_fast_zoomify_build"
    with h5py.File(path, "r+") as f:
        if temporary_group in f:
            del f[temporary_group]
        assembly = f[f"resolutions/{source}"].attrs.get("genome-assembly")
    new_offsets = initialize_level(path, temporary_group, names, lengths, target, assembly, args)
    old_offsets = np.r_[0, np.cumsum((lengths + source - 1) // source)]
    factor, nbins = target // source, int(new_offsets[-1])
    with tempfile.TemporaryDirectory(prefix="fast-zoomify-", dir=args.temp_dir) as scratch, h5py.File(path, "r+") as f:
        old, new = f[f"resolutions/{source}"], f[temporary_group]
        pixels = old["pixels"]
        row_index = old["indexes/bin1_offset"]
        nnz = count_sum = 0
        previous = -1
        last_log = time.monotonic()

        def append(records):
            nonlocal nnz, count_sum, previous
            if not len(records):
                return
            if int(records["key"][0]) <= previous or np.any(records["key"][1:] <= records["key"][:-1]):
                raise RuntimeError("Coarsened pixels are not globally sorted and unique")
            previous = int(records["key"][-1])
            count_sum += int(records["count"].sum())
            values = (records["key"] // np.uint64(nbins), records["key"] % np.uint64(nbins), records["count"])
            for field, data in zip(("bin1_id", "bin2_id", "count"), values):
                ds = new[f"pixels/{field}"]
                ds.resize((nnz + len(records),))
                ds[nnz:nnz + len(records)] = data
            nnz += len(records)

        for cid in range(len(names)):
            first, stop = map(int, old_offsets[cid:cid + 2])
            lo, end = int(row_index[first]), int(row_index[stop])
            while lo < end:
                probe = min(lo + args.chunksize, end)
                hi = end
                if probe < end:
                    old_row = int(pixels["bin1_id"][probe])
                    coarse_start = first + ((old_row - first) // factor) * factor
                    hi = int(row_index[coarse_start])
                    if hi <= lo:  # One coarse row exceeds the pixel budget: spill.
                        hi = int(row_index[min(coarse_start + factor, stop)])
                levels, serial = [], 0
                for a in range(lo, hi, args.chunksize):
                    b = min(a + args.chunksize, hi)
                    i, j = pixels["bin1_id"][a:b], pixels["bin2_id"][a:b]
                    ci = new_offsets[cid] + (i - first) // factor
                    cj_chrom = np.searchsorted(old_offsets, j, side="right") - 1
                    cj = new_offsets[cj_chrom] + (j - old_offsets[cj_chrom]) // factor
                    records = np.empty(b - a, dtype=RECORD)
                    records["key"] = ci.astype(np.uint64) * np.uint64(nbins) + cj.astype(np.uint64)
                    counts = pixels["count"][a:b]
                    if np.any(counts < 0) or (counts.dtype.kind == "u" and np.any(counts > np.iinfo(np.int64).max)):
                        raise ValueError("Counts must be nonnegative int64-compatible integers")
                    records["count"] = counts
                    records = aggregate(records)
                    if hi - lo <= args.chunksize:
                        append(records)
                    else:
                        run = os.path.join(scratch, f"{serial}.bin"); serial += 1
                        records.tofile(run)
                        level = 0
                        while level < len(levels) and levels[level] is not None:
                            dest = os.path.join(scratch, f"{serial}.bin"); serial += 1
                            merge_runs(levels[level], run, dest, args.io_chunk)
                            levels[level] = None
                            run = dest
                            level += 1
                        if level == len(levels):
                            levels.append(run)
                        else:
                            levels[level] = run
                runs = [x for x in levels if x is not None]
                if runs:
                    run = runs[0]
                    for other in runs[1:]:
                        dest = os.path.join(scratch, f"{serial}.bin"); serial += 1
                        merge_runs(run, other, dest, args.io_chunk)
                        run = dest
                    with open(run, "rb") as stream:
                        while True:
                            records = np.fromfile(stream, dtype=RECORD, count=args.io_chunk)
                            if not len(records):
                                break
                            append(records)
                    os.unlink(run)
                lo = hi
                if time.monotonic() - last_log > 10:
                    LOG.info("Coarsening %d -> %d: %s output pixels", source, target, f"{nnz:,}")
                    last_log = time.monotonic()
        if "sum" in old.attrs and count_sum != int(old.attrs["sum"]):
            raise RuntimeError("Coarsening failed count-sum conservation check")
        make_row_index(new, args.io_chunk)
        new.attrs.update({"nnz": nnz, "sum": count_sum, "fast-raw-complete": True})
        f.move(temporary_group, f"/resolutions/{target}")


def sparse_memory_estimate(nbins, nnz, chunksize):
    # Conservative planning estimate for CSR conversion, slicing, symmetric
    # expansion, vectors and allocator overlap; not an enforced OS memory cap.
    return 80 * nnz + 128 * nbins + 128 * chunksize


def load_upper(g, ignore_diags, chunksize):
    n, z = int(g.attrs["nbins"]), int(g.attrs["nnz"])
    index_dtype = np.int32 if max(n, z) <= np.iinfo(np.int32).max else np.int64
    pointers = g["indexes/bin1_offset"][:].astype(index_dtype, copy=False)
    indices = np.empty(z, dtype=index_dtype)
    values = np.empty(z, dtype=np.float64)
    for a in range(0, z, chunksize):
        b = min(a + chunksize, z)
        i = g["pixels/bin1_id"][a:b]
        j = g["pixels/bin2_id"][a:b]
        counts = g["pixels/count"][a:b]
        if np.any(i > j) or np.any(j >= n) or np.any(i < 0) or np.any(counts < 0):
            raise ValueError("Invalid symmetric-upper count pixels")
        indices[a:b] = j
        values[a:b] = counts
        if ignore_diags:
            values[a:b][j - i < ignore_diags] = 0
    upper = csr_matrix((values, indices, pointers), shape=(n, n), copy=False)
    upper.eliminate_zeros()
    return upper


def blacklist_mask(path, names, lengths, offsets, resolution, keep):
    if not path:
        return
    ids = {name: cid for cid, name in enumerate(names)}
    with open(path) as stream:
        for line in stream:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError("Blacklist must be a BED file with at least three columns")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in ids:
                continue
            cid = ids[chrom]
            start, end = max(0, start), min(int(lengths[cid]), end)
            if start < end:
                a = int(offsets[cid]) + start // resolution
                b = int(offsets[cid]) + (end + resolution - 1) // resolution
                keep[a:b] = False


def initial_keep(upper, names, lengths, offsets, resolution, args):
    n = upper.shape[0]
    keep = np.ones(n, dtype=bool)
    if args.min_nnz > 0:
        # Counts both ends of each retained pixel, including twice on diagonal,
        # matching Cooler's marginalization when --ignore-diags=0.
        degree = np.diff(upper.indptr).astype(np.int64)
        degree += np.bincount(upper.indices, minlength=n)
        keep[degree < args.min_nnz] = False
        del degree
    ones = np.ones(n)
    marg = upper @ ones
    marg += upper.T @ ones
    del ones
    if args.min_count:
        keep[marg < args.min_count] = False
    if args.mad_max > 0:
        for lo, hi in zip(offsets[:-1], offsets[1:]):
            part = marg[int(lo):int(hi)]
            positive = part[part > 0]
            if len(positive):
                part /= np.median(positive)
            else:
                # Reproduce Cooler's empty-chromosome normalization semantics.
                part[:] = np.nan
        positive = marg[marg > 0]
        if len(positive):
            logged = np.log(positive)
            median = np.median(logged)
            deviation = np.median(np.abs(logged - median))
            cutoff = np.exp(median - args.mad_max * deviation)
            keep[marg < cutoff] = False
    blacklist_mask(args.blacklist, names, lengths, offsets, resolution, keep)
    return keep


class ParallelMatvec:
    """CSR row views share matrix buffers; threads return only their row range."""
    def __init__(self, matrix, threads):
        self.matrix = matrix
        self.pool = None
        self.parts = []
        n = matrix.shape[0]
        if threads > 1 and n and matrix.nnz:
            targets = np.linspace(0, matrix.nnz, threads + 1, dtype=np.int64)
            edges = np.unique(np.r_[0, np.searchsorted(matrix.indptr, targets), n])
            edges = edges[(edges >= 0) & (edges <= n)]
            for lo, hi in zip(edges[:-1], edges[1:]):
                a, b = matrix.indptr[lo], matrix.indptr[hi]
                part = csr_matrix((int(hi - lo), n), dtype=np.float64)
                # Assign public CSR buffers directly to prevent SciPy from
                # copying/downcasting large index slices in the constructor.
                part.data = matrix.data[a:b]
                part.indices = matrix.indices[a:b]
                part.indptr = matrix.indptr[lo:hi + 1] - a
                part.has_sorted_indices = True
                self.parts.append((int(lo), int(hi), part))
            self.pool = ThreadPoolExecutor(max_workers=threads)

    def __call__(self, weights):
        if self.pool is None:
            return self.matrix @ weights
        result = np.empty(len(weights), dtype=np.float64)
        futures = [(lo, hi, self.pool.submit(part.dot, weights)) for lo, hi, part in self.parts]
        for lo, hi, future in futures:
            result[lo:hi] = future.result()
        return result

    def close(self):
        if self.pool is not None:
            self.pool.shutdown()


def ice(matrix, args):
    bias = np.ones(matrix.shape[0], dtype=np.float64)
    scale, variance, converged, reason = np.nan, np.nan, False, "iteration limit"
    iterations = 0
    multiply = ParallelMatvec(matrix, args.threads)
    try:
        for iteration in range(args.max_iters):
            iterations = iteration + 1
            marg = multiply(bias)
            marg *= bias
            nz = marg[marg != 0]
            if not len(nz):
                bias[:] = np.nan
                variance, reason = np.nan, "no usable contacts after filtering"
                break
            if not np.all(np.isfinite(nz)):
                reason = "non-finite marginals"
                break
            scale, variance = float(nz.mean()), float(nz.var())
            if not math.isfinite(scale) or scale <= 0 or not math.isfinite(variance):
                reason = "non-finite or nonpositive marginal statistics"
                break
            marg /= scale
            marg[marg == 0] = 1
            bias /= marg
            if not np.all(np.isfinite(bias)):
                reason = "non-finite weights"
                break
            if iteration == 0 or iterations % 10 == 0 or variance < args.tol:
                LOG.info("ICE iteration %d: variance=%.8g", iterations, variance)
            if variance < args.tol:
                converged, reason = True, "converged"
                break
    finally:
        multiply.close()
    if math.isfinite(scale) and scale > 0:
        bias /= np.sqrt(scale)
    else:
        bias[:] = np.nan
    return bias, dict(scale=scale, var=variance, converged=converged,
                      iterations=iterations, reason=reason)


def balance_level(path, resolution, names, lengths, args):
    group = f"resolutions/{resolution}"
    with h5py.File(path, "r") as f:
        g = f[group]
        n, z = int(g.attrs["nbins"]), int(g.attrs["nnz"])
        estimate = sparse_memory_estimate(n, z, args.chunksize)
        LOG.info("%d bp ICE: %s bins, %s pixels; planning estimate %.1f GiB",
                 resolution, f"{n:,}", f"{z:,}", estimate / GIB)
        if estimate > args.memory_gb * GIB:
            raise MemoryError(f"Planning estimate {estimate/GIB:.1f} GiB exceeds --memory-gb {args.memory_gb:g}")
        offsets = g["indexes/chrom_offset"][:]
        ignore_diags = max(args.ignore_diags, math.ceil(args.ignore_dist / resolution))
        upper = load_upper(g, ignore_diags, args.chunksize)
    keep = initial_keep(upper, names, lengths, offsets, resolution, args)
    retained = int(keep.sum())
    LOG.info("%d bp: retained %s / %s bins (%.2f%%)", resolution, f"{retained:,}", f"{n:,}", 100*retained/n)
    # Compact only the solve; raw pixels and the complete bin grid are retained.
    if retained < n:
        upper = upper[keep, :][:, keep]
    # Full symmetric CSR allows row-parallel matvecs with no scatter races.
    # U + U.T deliberately doubles diagonal contributions, like Cooler ICE.
    matrix = upper + upper.T
    del upper
    gc.collect()
    bias, stats = ice(matrix, args)
    del matrix
    stats.update(tol=args.tol, min_nnz=args.min_nnz, min_count=args.min_count,
                 mad_max=args.mad_max, cis_only=False, ignore_diags=ignore_diags,
                 divisive_weights=False, retained_bins=retained,
                 algorithm=f"fast_zoomify-{VERSION}")
    if not stats["converged"]:
        LOG.warning("%d bp ICE did not converge: %s", resolution, stats["reason"])
        if args.convergence_policy == "error":
            raise RuntimeError(f"ICE did not converge at {resolution} bp; output is resumable")
        if args.convergence_policy == "store_nan":
            bias[:] = np.nan
        elif stats["reason"] != "iteration limit" or not np.all(np.isfinite(bias)):
            LOG.warning("%d bp: no finite final ICE result to retain; writing NaN weights", resolution)
            bias[:] = np.nan
    with h5py.File(path, "r+") as f:
        bins = f[group + "/bins"]
        temporary = "_fast_weight"
        if temporary in bins:
            del bins[temporary]
        ds = bins.create_dataset(temporary, shape=(n,), dtype="f8",
                                 chunks=(min(n, args.io_chunk),), compression="gzip", compression_opts=args.compression)
        cursor = 0
        finite = 0
        for a in range(0, n, args.io_chunk):
            b = min(a + args.io_chunk, n)
            mask = keep[a:b]
            count = int(mask.sum())
            values = np.full(b - a, np.nan)
            values[mask] = bias[cursor:cursor + count]
            cursor += count
            finite += int(np.isfinite(values).sum())
            ds[a:b] = values
        stats["finite_weights"] = finite
        ds.attrs.update(stats)
        if "weight" in bins:
            del bins["weight"]
        bins.move(temporary, "weight")
    del bias, keep
    gc.collect()
    return stats


def file_fingerprint(path):
    p = Path(path).resolve()
    stat = p.stat()
    return dict(path=str(p), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def prepare_output(args, source, resolutions):
    output = str(Path(args.out).resolve())
    if output == str(Path(source.filename).resolve()):
        raise ValueError("Output must be a different file from input")
    manifest = dict(version=VERSION, source=file_fingerprint(source.filename), group=source.root,
                    resolutions=resolutions, ignore_diags=args.ignore_diags, ignore_dist=args.ignore_dist,
                    min_nnz=args.min_nnz, min_count=args.min_count, mad_max=args.mad_max, tol=args.tol,
                    blacklist=file_fingerprint(args.blacklist) if args.blacklist else None)
    encoded = json.dumps(manifest, sort_keys=True)
    if os.path.exists(output):
        if not args.resume:
            raise FileExistsError("Output exists; use --resume with matching input/resolutions/filter settings")
        with h5py.File(output, "r+") as f:
            if f.attrs.get("fast-zoomify-manifest") != encoded:
                raise ValueError("Resume settings or source fingerprint differ; choose a new output")
            if "_fast_zoomify_build" in f:
                del f["_fast_zoomify_build"]
        return output
    with h5py.File(output, "x") as f:
        f.attrs.update({"format": "HDF5::MCOOL", "format-version": 2,
                        "fast-zoomify-manifest": encoded, "fast-zoomify-status": "incomplete"})
        f.create_group("resolutions")
    return output


def run(args):
    source, names, lengths, offsets = source_info(args.input)
    resolutions = resolutions_for(args.resolutions, source.binsize, int(lengths.sum()))
    if args.temp_dir:
        Path(args.temp_dir).mkdir(parents=True, exist_ok=True)
    path = prepare_output(args, source, resolutions)
    LOG.info("Resolutions (all will be balanced): %s", ",".join(map(str, resolutions)))
    try:
        with h5py.File(path, "r+") as out:
            base_group = f"resolutions/{source.binsize}"
            if base_group not in out:
                LOG.info("Copying the raw %d bp base matrix", source.binsize)
                if "_fast_base_copy" in out:
                    del out["_fast_base_copy"]
                with source.open() as src:
                    src.file.copy(src, out, name="_fast_base_copy")
                if "weight" in out["_fast_base_copy/bins"]:
                    del out["_fast_base_copy/bins/weight"]
                out["_fast_base_copy"].attrs["fast-raw-complete"] = True
                out.move("_fast_base_copy", base_group)
        # Generate all raw levels before attempting potentially difficult ICE.
        for resolution in resolutions[1:]:
            with h5py.File(path, "r") as f:
                if f"resolutions/{resolution}" in f:
                    continue
            predecessor = max(r for r in resolutions if r < resolution and resolution % r == 0)
            LOG.info("Coarsening %d -> %d", predecessor, resolution)
            coarsen_level(path, predecessor, resolution, names, lengths, args)
        failed = []
        for resolution in resolutions:
            with h5py.File(path, "r") as f:
                key = f"resolutions/{resolution}/bins/weight"
                if key in f and bool(f[key].attrs.get("converged", False)):
                    LOG.info("Skipping already converged %d bp level", resolution)
                    continue
            stats = balance_level(path, resolution, names, lengths, args)
            if not stats["converged"]:
                failed.append(resolution)
        with h5py.File(path, "r+") as f:
            f.attrs["fast-zoomify-status"] = "complete" if not failed else "nonconverged"
            f.attrs["fast-zoomify-failed-resolutions"] = json.dumps(failed)
            if "fast-zoomify-error" in f.attrs:
                del f.attrs["fast-zoomify-error"]
        if failed:
            LOG.error("Nonconverged levels: %s. Inspect weight attributes; --resume can retry with more iterations.", failed)
            return 2
        LOG.info("Complete: %s; all %d levels converged", path, len(resolutions))
        return 0
    except Exception as exc:
        with h5py.File(path, "r+") as f:
            f.attrs["fast-zoomify-status"] = "incomplete"
            f.attrs["fast-zoomify-error"] = str(exc)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="Input fixed-bin .cool path or Cooler URI")
    parser.add_argument("--out", "-o", required=True, help="New output .mcool filename")
    parser.add_argument("--resolutions", "-r", default="B", help="B: binary; N: 1-2-5 progression; or comma-separated bp sizes. Base always included")
    parser.add_argument("--threads", type=int, default=8, help="Shared-memory ICE matvec threads; coarsening uses one process")
    parser.add_argument("--chunksize", type=int, default=2_000_000, help="Input pixels per coarsening/loading batch")
    parser.add_argument("--io-chunk", type=int, default=250_000, help="Output/index/spill-merge block length")
    parser.add_argument("--compression", type=int, choices=range(10), default=1)
    parser.add_argument("--memory-gb", type=float, default=128, help="GiB planning budget for ICE; conservative estimate, not a hard OS limit")
    parser.add_argument("--temp-dir", help="Scratch for unusually dense coarse rows; default system temp")
    parser.add_argument("--ignore-diags", type=int, default=2)
    parser.add_argument("--ignore-dist", type=int, default=0, help="Ignored distance in bp; combined with --ignore-diags as in Cooler")
    parser.add_argument("--min-nnz", type=int, default=10)
    parser.add_argument("--min-count", type=int, default=0)
    parser.add_argument("--mad-max", type=int, default=5)
    parser.add_argument("--blacklist", help="Headerless BED; overlapping bins excluded independently at each resolution")
    parser.add_argument("--tol", type=float, default=1e-5)
    parser.add_argument("--max-iters", type=int, default=200)
    parser.add_argument("--convergence-policy", choices=("store_nan", "error", "store_final"), default="store_nan",
                        help="Retain NaNs (default), retain finite final weights at the iteration limit, or abort without writing weights.")
    parser.add_argument("--resume", action="store_true", help="Reuse this script's output with the same source/resolutions/filters; retry failures")
    args = parser.parse_args(argv)
    for name in ("threads", "chunksize", "io_chunk", "max_iters", "memory_gb", "tol"):
        if not math.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    for name in ("ignore_diags", "ignore_dist", "min_nnz", "min_count", "mad_max"):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be nonnegative")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        return run(args)
    except (ValueError, OSError, RuntimeError, MemoryError) as exc:
        LOG.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
