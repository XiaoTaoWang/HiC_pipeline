# Efficient fixed-bin zoomification and genome-wide ICE

`fast_zoomify.py` is a standalone command. It imports Cooler and its existing
NumPy, pandas, SciPy, and h5py dependencies; it does not change the Cooler package
and does not require the earlier pairs-loading script or pypairix.

It generates a standard `.mcool` and balances **every requested resolution,
including the original 10 bp level**, independently from raw counts.

## Recommended command for your machine

For 64 available CPUs and 256 GB RAM:

```bash
python fast_zoomify.py out.10bp.cool \
  --out out.mcool \
  --resolutions B \
  --threads 8 \
  --memory-gb 160 \
  --chunksize 2000000 \
  --io-chunk 250000 \
  --compression 1 \
  --max-iters 200 \
  --convergence-policy store_nan \
  --temp-dir /path/to/scratch
```

Replace the scratch path. Prefer node-local SSD when available. Scratch is used
only when one coarse matrix row exceeds the pixel batch size. Most coarsening
batches are aggregated and written without temporary files. For shared-storage
input/output, staging the input and output on local disk can also help; copy the
completed output back afterward.

`B` matches Cooler's default binary progression: 10, 20, 40, 80, ... up to the
largest level at or below approximately genome length / 256. For primary hg38,
the final binary level is 10,485,760 bp. This produces 21 resolutions including
10 bp. Choose only scientifically useful levels if a smaller pyramid suffices.
All levels in the chosen pyramid will still be balanced.

Alternatives:

```bash
# 10, 20, 50, 100, 200, 500, 1000, ...
python fast_zoomify.py out.10bp.cool -o out.nice.mcool --resolutions N --memory-gb 160

# Explicit levels; the original 10 bp level is always included and balanced.
python fast_zoomify.py out.10bp.cool -o out.selected.mcool \
  --resolutions 10,50,100,500,1000,5000,10000 --memory-gb 160
```

The script requires Python 3.10 or newer. Run it in your Cooler environment.
The output parent directory must exist. Existing output is never overwritten;
use `--resume` for output previously created by this script.

## Algorithm and compatibility

### Coarsening

- Reads sorted pixels in bounded batches and computes coarser IDs directly from
  chromosome bin offsets; no genome-wide coordinate join or source bin table.
- Respects chromosome boundaries and partial final bins.
- Aligns batches to coarse rows. Oversized rows use sorted spill runs and bounded,
  vectorized merging, so duplicate coarse pixels are combined across batches.
- Uses the closest available finer resolution that divides the target exactly.
- Uses `cooler.create_cooler()` for small schema templates, then streams large
  bin tables, pixels, and row indexes using h5py.
- Preserves int64 raw counts and checks count-sum conservation when the input
  provides a `sum` attribute. Coarsened counts are never derived from balanced data.

Coarsening currently uses one process. `--threads` controls ICE sparse
matrix–vector operations, not concurrent HDF5 writers. The vectorized coarsener
removes the principal bin-table and pandas-grouping overhead; its actual
throughput depends on storage and the number of requested levels.

### ICE

- Reads the upper-triangular pixels once for a resolution and applies diagonal
  exclusions.
- Computes the same initial genome-wide minimum-nonzero, minimum-count, MAD-max,
  and optional blacklist filters as Cooler 0.10.4.
- Compacts the retained bin IDs for the solve only; the output still contains
  the complete bin grid and all raw pixels.
- Builds a symmetric sparse CSR matrix in RAM. Threads share its arrays and the
  current weight vector. Each thread returns only its own matrix-row range.
- Performs the Cooler-style multiplicative ICE update, raw marginal-variance
  convergence check, and final scale normalization.
- Expands weights back to the original bin grid, with NaNs for excluded bins,
  and writes standard `bins/weight` datasets with convergence attributes.

The full bin table is not read every iteration, and genome-length marginal
vectors are not serialized for every pixel chunk. No dense square matrix is
constructed. There is deliberately only one active resolution in memory.

This is a numerically compatible implementation of genome-wide ICE, not a call
to Cooler's original out-of-core balancing loop. Floating-point summation order
differs, so bitwise-identical weights are not promised. Tests compare weights and
scale/variance statistics at tight numerical tolerances. With ignored diagonals
set to zero, diagonal contributions are counted twice during ICE, matching
Cooler's current marginalization convention.

## Tuning and memory

| Option | Starting value | Meaning |
|---|---:|---|
| `--threads` | 8 | Shared-memory CSR matvec workers. Compare 4 and 8; 64 is unlikely to help a memory-bandwidth-bound solve. |
| `--memory-gb` | 160 | Planning budget in **GiB**, checked before each ICE matrix allocation. Not a hard memory cap or preallocation. |
| `--chunksize` | 2,000,000 | Pixels per input/coarsening batch. Larger values can reduce spill work but increase temporary sorting memory. |
| `--io-chunk` | 250,000 | HDF5, row-index, and spill-merge batch length. |
| `--compression` | 1 | Gzip level for newly generated levels and weights. The copied base retains its original compression. |
| `--max-iters` | 200 | Iteration limit per resolution. More iterations help only when a solve is approaching convergence. |
| `--tol` | 1e-5 | Same marginal-variance criterion as Cooler. Changing it changes the numerical requirement. |

For approximately 309 million bins and one billion pixels, the script's
conservative allocation estimate is approximately **112 GiB**, before knowing
how much filtering will reduce the matrix. Thus `--memory-gb 160` is a reasonable
initial setting on your 256 GB node. This is an engineering estimate, not a
measured full-scale peak or an enforced OS allocation limit. Request enough
memory from your scheduler as well; this flag does not reserve cluster RAM.

Reducing `--chunksize` cannot make the CSR matrix itself smaller. If the planning
guard fails, the output's completed raw levels remain available. A genuinely
out-of-core balancer would be needed if the sparse matrix and workspace cannot
fit. Simply increasing the flag is appropriate only when sufficient physical
and scheduler-allocated memory is available.

There is no validated wall-time estimate for the full billion-contact pyramid.
Time depends on pixel reduction across levels, convergence iterations, retained
coverage, memory bandwidth, and shared-storage throughput. Many fine levels can
retain almost as many pixels as the base, increasing both output size and work.

## Scientific settings and convergence

Defaults match Cooler genome-wide ICE:

```text
--ignore-diags 2
--ignore-dist 0
--min-nnz 10
--min-count 0
--mad-max 5
--tol 1e-5
--max-iters 200
```

`--ignore-dist` is in bp; each resolution uses the greater of `--ignore-diags`
and `ceil(ignore_dist / resolution)`. `--blacklist` accepts a headerless BED file,
and bins overlapping its intervals are excluded independently at each level.
Do not lower filtering thresholds solely to manufacture finite weights at 10 bp.

Balancing every level means that every level is **attempted**. Sparse or poorly
supported levels may have many excluded bins or fail to converge. Additional
CPUs and iterations cannot guarantee usable ICE weights.

The default failure policy is intentionally stricter than Cooler's default:

- `store_nan` (default): write all-NaN weights for a failed level, continue other
  levels, then exit with status **2** if any level failed.
- `error`: stop at the first failed level with exit status **1**; completed raw
  levels and previously converged weights remain available.
- `store_final`: store the last finite weights at the iteration limit, mark `converged=False`, and
  continue; final exit status is still **2** if a level failed. These weights
  must not be treated as converged results. If ICE encounters invalid numerics
  rather than reaching the iteration limit, this policy writes NaN weights.

Exit status **0** means every requested level converged. Other input, memory,
or I/O errors produce status **1**. If no usable contacts survive filtering,
this script marks the level nonconverged; this deliberately avoids Cooler's
empty-marginal edge case reporting success with all-NaN weights.

Inspect convergence without loading large weight vectors:

```python
import h5py

with h5py.File("out.mcool", "r") as f:
    print("Job status:", f.attrs.get("fast-zoomify-status"))
    for resolution in sorted(f["resolutions"], key=int):
        bins = f[f"resolutions/{resolution}/bins"]
        if "weight" not in bins:
            print(resolution, "not yet balanced")
            continue
        a = bins["weight"].attrs
        print(resolution, "converged:", a["converged"],
              "iterations:", a["iterations"], "finite weights:", a["finite_weights"],
              "reason:", a["reason"])
```

Use normal Cooler access afterward:

```python
import cooler

c = cooler.Cooler("out.mcool::resolutions/10")
raw = c.matrix(balance=False, sparse=True).fetch("chr1:1000000-1010000")
# Check convergence before interpreting balanced values.
balanced = c.matrix(balance=True, sparse=True).fetch("chr1:1000000-1010000")
```

## Resume

All raw levels are generated before balancing starts. The output file remains
available after a controlled failure; a root status attribute records completion.
To retry a slow-converging level with more iterations:

```bash
python fast_zoomify.py out.10bp.cool \
  --out out.mcool --resolutions B \
  --threads 8 --memory-gb 160 \
  --max-iters 500 --convergence-policy store_nan \
  --temp-dir /path/to/scratch --resume
```

Resume skips completed raw levels and already converged weight datasets. Failed
levels restart ICE from unit weights; this is not an iteration checkpoint.
You can change thread count, batch sizes, memory budget, iteration limit, or
convergence policy. Input fingerprint, resolutions, tolerance, blacklist, and
filter settings must match. Choose a new output when changing scientific settings.

Do not run two writers against the same output at the same time.

## Scope and validation

Supported: complete fixed-bin grids, symmetric-upper storage, nonnegative integer
counts, chromosomes no longer than 2**31 - 1 bp, at most 2**32 bins, and count sums
that fit int64. This includes the output of `fast_pairix_to_cool.py` for primary
hg38. Variable-width bins, square storage, floating-point count matrices, cis-only
balancing, and trans-only balancing are not implemented. Input is assumed to be
a valid Cooler; this script is not a full corruption audit.

Tested on Python 3.12.14, Cooler 0.10.4, NumPy 2.3.5, SciPy 1.18.1, pandas 2.2.3,
and h5py 3.16.0 on macOS ARM64. Linux and other dependency versions have not been
tested here. No newly installed dependencies are needed beyond Cooler's dependencies.

Nine automated tests passed, covering:

- Exact coarsened pixels against Cooler, including forced multi-run spills,
  chromosome boundaries, partial final bins, and explicit zero-count pixels.
- ICE weights matching `cooler.balance_cooler()` at relative tolerance 1e-10,
  with one/multiple threads, diagonal inclusion/exclusion, coverage filters,
  MAD filtering, and blacklist intervals.
- Balanced regional matrix queries, empty levels, iteration-limit failures,
  resume/retry, incompatible-resume rejection, memory guard, and uneven thread
  partitions.

A convergent synthetic **10 bp** benchmark with 100,000 bins and approximately
two million pixels used one compute thread/process in both implementations:

| Method | Elapsed | Peak process RSS |
|---|---:|---:|
| This script, base level only (including copy and ICE) | 0.50 s | 222 MB |
| Stock `cooler balance`, one worker | 1.23 s | 145 MB |

99,741 bins received finite weights, ICE converged in seven iterations, and
weights matched within relative tolerance 1e-10. The optimized implementation
trades RAM for fewer reads and more efficient iteration. These are single,
synthetic, local-storage measurements, not a forecast for one billion pixels.

A separate five-million-contact, ten-million-bin test generated 1 kb and 5 kb
levels and balanced those successfully in a total pipeline run of 2.81 seconds
with 806 MB peak RSS. Its 10 bp level had no usable contacts after filtering and
was correctly flagged as failed. That run is evidence for failure handling and
streaming coarsening, **not** a successful large 10 bp ICE benchmark.
