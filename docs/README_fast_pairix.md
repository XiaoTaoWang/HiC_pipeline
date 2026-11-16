# Standalone high-resolution Pairix → Cooler loader

`fast_pairix_to_cool.py` reads an existing sorted, BGZIP-compressed, Pairix-indexed
4DN pairs file and writes a normal, unbalanced, symmetric-upper `.cool` file.
It does not modify your Cooler installation.

## Run

Use the Python environment in which Cooler and pypairix are installed. For your
64 available CPUs, 256 GB RAM, and approximately 10 GB compressed input on shared
storage, start with this configuration:

```bash
python fast_pairix_to_cool.py \
  --assembly hg38 \
  --nproc 8 \
  --window-bp 5000000 \
  --buffer-contacts 4000000 \
  --io-chunk 250000 \
  --compression 1 \
  --temp-dir /path/to/scratch \
  hg38.chrom.sizes:10 in.pairs.gz out.10bp.cool
```

Replace the scratch path with node-local SSD scratch if available; otherwise,
use the fastest available shared scratch directory. Omitting
`--temp-dir` uses the system temporary directory. The output's parent directory
must already exist. The input `.px2` index must already exist and match the pairs
file. Install `pypairix` separately if your Cooler environment does not include it.
The script requires Python 3.10 or newer, Cooler, pypairix, numpy, pandas, and h5py;
the last three are Cooler dependencies. No pysam dependency is needed to run it.

For subsequent use:

```bash
cooler info out.10bp.cool
cooler dump -t pixels out.10bp.cool | head
```

```python
import cooler

c = cooler.Cooler("out.10bp.cool")
print(c.info)
# Keep matrix requests regional at this resolution.
region = c.matrix(balance=False, sparse=True).fetch("chr1:1000000-1010000")
```

## What is different

1. **Interval queries:** each worker queries a large, bin-aligned interval against
   all relevant partner chromosomes, rather than querying every 10 bp bin.
2. **Numeric buffers:** contacts become uint64 pixel keys. Sorting and counting
   happen in NumPy, with no DataFrame per occupied row.
3. **Bounded aggregation:** filled buffers are written as sorted numeric runs.
   Dense intervals trigger bounded, vectorized merges. Even one extremely dense
   bin does not require retaining all its contacts in RAM.
4. **Small worker state:** workers receive chromosome lengths, offsets, options,
   and the input path. They do not receive the full bin table. Each worker opens
   its own Pairix handle once.
5. **Bounded scheduling:** at most twice the worker count of intervals are
   submitted. Workers return paths to numeric runs rather than large arrays.
6. **One output writer:** globally sorted, unique pixels are appended by a single
   HDF5 writer. Bins and row offsets, including empty rows, are written in blocks.

The script uses the public `cooler.create_cooler()` API to create a small schema
template, then replaces its large datasets using h5py. It uses `cooler.Cooler`
and `cooler.fileops.is_cooler()` for final structural checks. This direct streaming
step is necessary to avoid the full-genome pandas bin table required by Cooler's
standard creation interface. It does not use private Cooler modules.

## Tuning

| Option | Default | Effect |
|---|---:|---|
| `--nproc` | 4 | Worker count; raises aggregate worker/index memory |
| `--window-bp` | 1,000,000 | Genomic interval width, rounded down to a bin multiple |
| `--buffer-contacts` | 2,000,000 | Maximum input contacts in a worker's numeric buffer |
| `--io-chunk` | 250,000 | Merge, index, and HDF5 batch length |
| `--compression` | 1 | Gzip level, from 0 through 9 |
| `--zero-based` | off | Use only when both coordinates and index are zero-based |
| `--force` | off | Replace existing output only after successful completion |

The table lists the script's built-in defaults; the command above overrides them
for your machine. Start with that command, compare 8 versus 16 workers, then
compare 1 Mb versus 5 Mb intervals using completed contacts per second. Leave
`--io-chunk 250000` and `--compression 1` unchanged initially. If RAM is tight,
lower `--buffer-contacts` to 500000
and/or lower the worker count. This increases spill/merge work. Larger windows
reduce indexed queries but can increase temporary I/O and task imbalance. If
many intervals spill, reducing `--window-bp` can help. The buffer cap still
protects memory for uneven coverage. Unlike `--max-split`, buffer size is a
direct limit on the number of contacts retained before aggregation.

Memory scales mainly with the number of workers, numeric buffer size, merge
buffers, and per-worker Pairix index, rather than the full number of genomic
bins. The recommended four-million-contact key buffer is 32 MB per worker
(the built-in default is 16 MB), but sorting, unique/count
arrays, and Python/library overhead make actual process memory considerably
larger. Do not interpret the raw buffer size as a process memory limit.

Prefer node-local disk for scratch when available. Each numeric run uses 16 bytes per distinct pixel
in that run. Spill merging temporarily holds input and output runs together;
duplicates can occur in several intermediate runs before merging. Completed
interval files are removed as the writer consumes them. For a billion contacts,
allow tens of GB of scratch headroom, plus space for the output. Very dense
intervals and slow or network-mounted scratch can substantially affect runtime.

There is still O(number of bins) output work: a conventional full-genome 10 bp
Cooler must store every bin and every row offset. It is streamed rather than
held in memory. Downstream Cooler commands have their own memory requirements.

## Input and output contract

- Fixed-width bins supplied as `chrom.sizes:BIN_SIZE`; BED/variable bins are not
  supported by this script.
- Point pairs with one record per observed contact, using a consistent triangular
  representation. Cis pairs must satisfy `pos1 <= pos2`. Contacts are counted,
  not deduplicated; PCR/read duplicates must be removed upstream if desired.
- Standard 4DN one-based coordinates are the default. The index and coordinates
  must agree. Zero-based input has been separately tested.
- Input chromosome-block order may differ from chrom.sizes order, including
  reversed trans blocks. The output uses chrom.sizes order. Both orientations
  of the same trans chromosome pair are rejected.
- Contacts involving chromosomes absent from chrom.sizes are ignored, with a
  warning. Chromosomes present in chrom.sizes but absent from the input remain
  in the output as empty bins.
- Only indexed, in-range point pairs are read. This is not a whole-input audit
  for malformed or out-of-bounds records that are absent from queried regions.
- Partial final bins are retained. Interval boundaries are aligned to bins so
  that no output pixel is divided between intervals.
- Output `bin1_id`, `bin2_id`, and `count` are int64. Int64 counts are valid in
  Cooler and avoid the int32 count ceiling. Starts/ends are int32. Each chromosome
  must have length at most 2**31 - 1; the complete bin grid must have at most
  2**32 bins. Primary hg38 at 10 bp is within these limits.
- Output is a standalone `.cool`, not an appended group inside an existing file.
  The assembly label is retained. Other pairs-header metadata is not copied.
- A private partial output is published only after successful completion and
  validation. With `--force`, a previous output is preserved if the run fails.
- Sorting, indexing, mapping-quality filtering, deduplication, balancing, and
  multiresolution generation are not performed by this loader.

## Validation and measured examples

Tested with Python 3.12.14, Cooler 0.10.4, pypairix 0.3.9, NumPy 2.3.5,
pandas 2.2.3, and h5py 3.16.0 on macOS ARM64. The script uses portable spawned
workers; Linux and other package versions have not been exercised here.

Nine automated tests passed. They cover exact counts against the stock Pairix
loader and an independent count calculation; one/two workers; one/zero-based
coordinates; chromosome-order differences; bin and interval boundaries;
duplicates crossing many forced spills; empty chromosomes; empty output;
single-bin output; randomized spill merging; row indexes; matrix/region queries;
coarsening; invalid triangular input; and preservation of existing output on
failure. A separate ten-million-bin sparse output passed regional and long-gap
row-index checks.

Small synthetic benchmark, **100,000 contacts, 80,000 bins, 10 bp resolution**:

| Loader | Elapsed time | Peak process RSS |
|---|---:|---:|
| Stock `cooler cload pairix --nproc 1` | 49.9 s | 222 MB |
| This script, `--nproc 1 --compression 6` | 0.49 s | 122 MB |

Both used gzip level 6 for this comparison. Pixel IDs and counts matched exactly;
the count dtype differs (stock int32, script int64). Measurements include Python
imports, loading, output writing, and indexing. This is one synthetic run of
each implementation on local storage; it is not a billion-contact benchmark or
a promise of a 100-fold production speedup.

A separate **10,000,000-bin, two-contact** case completed in 0.66 seconds with
143 MB peak RSS using one worker and default compression. It demonstrates
streaming empty-bin metadata, not realistic contact-processing throughput.

Benchmark full-density representative intervals on your own CPU and storage
before setting the production wall-time or memory request. The billion-contact
case has not been run here.

### Larger tuning benchmark

A subsequent synthetic test used 5,000,000 pairs, three chromosomes totaling
100 Mb (10,000,000 bins at 10 bp), a 50.4 MB compressed input, and local scratch.
It included cis contacts with varied genomic distances and trans contacts.
Each configuration produced exactly the same pixel table (4,985,744 pixels).

| Workers | Interval width | Contact buffer | Elapsed | Peak summed process RSS |
|---:|---:|---:|---:|---:|
| 1 | 1 Mb | 2 million | 5.66 s | 179 MB |
| 4 | 1 Mb | 2 million | 2.89 s | 577 MB |
| 4 | 10 Mb | 2 million | 3.10 s | 729 MB |
| 4 | 10 Mb | 250,000 | 3.13 s | 827 MB |

RAM was sampled across the parent and worker processes, including shared pages
in each process's RSS. Measurements include startup and final indexing. The
files were on local storage and could be cached; these are single runs, not
repeated statistical benchmarks. The three-chromosome distribution is not a
full hg38 workload. The smaller contact buffer triggered additional merge
allocations, so observed peak RSS was not monotonically lower with buffer size.

For the user's proposed machine (64 available CPUs, 256 GB RAM, ~10 GB compressed
input on shared storage), an initial configuration to try is 8 workers, 5 Mb
intervals, a 4-million-contact buffer, the default 250,000-record I/O batches,
and gzip level 1. Compare 8 versus 16 workers and 1 versus 5 Mb intervals before
increasing further. Node-local scratch, when available, avoids sending the
intermediate runs over the shared filesystem.

A provisional production allowance for one billion contacts is 30 minutes to
3 hours with reasonably responsive shared storage, with 3–6+ hours possible
under contention or substantial indexed rereading/spilling. This is an engineering
estimate, not an extrapolated benchmark or a bound. Expect a few GB of process
RAM under ordinary coverage; reserving 16–32 GB for the first 8-worker run is
prudent. Confirm both on the actual input and filesystem.

## Source references

- [Cooler creation API](https://cooler.readthedocs.io/en/stable/api.html#cooler.create_cooler)
- [Cooler storage schema](https://cooler.readthedocs.io/en/stable/schema.html)
- [Original Pairix aggregator](https://cooler.readthedocs.io/en/stable/_modules/cooler/create/_ingest.html)
- [Pairix input and query documentation](https://github.com/4dn-dcic/pairix)
