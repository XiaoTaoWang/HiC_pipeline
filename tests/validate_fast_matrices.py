"""Real numerical smoke test. Requires the runHiC environment, bgzip and pairix.

Run from the source tree: python tests/validate_fast_matrices.py --work-dir /scratch/new-dir
The work directory must not exist; outputs are kept for inspection.
"""
import argparse
from collections import Counter
from pathlib import Path
import subprocess
import sys

import cooler
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from runHiC.binning import mcool_from_pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', required=True)
    args = parser.parse_args()
    work = Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=False)
    chroms = [('chr1', 103), ('chr2', 87)]
    bins = [(name, start + 1) for name, length in chroms for start in range(0, length, 10)]
    rows = []
    for i, (a, x) in enumerate(bins):
        for j in range(i, len(bins)):
            b, y = bins[j]
            for repeat in range(1 + (i * 7 + j * 3) % 4):
                rows.append((a, b, x, y, f'r{i}_{j}_{repeat}'))
    rows.sort()
    pairs_text = work / 'toy.pairs'
    with pairs_text.open('w', encoding='utf-8', newline='\n') as out:
        out.write('## pairs format v1.0.0\n#shape: upper triangle\n#genome_assembly: toy\n')
        for name, length in chroms:
            out.write(f'#chromsize: {name} {length}\n')
        out.write('#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n')
        for a, b, x, y, read_id in rows:
            out.write(f'{read_id}\t{a}\t{x}\t{b}\t{y}\t+\t-\n')
    pairs = work / 'toy.pairs.gz'
    with pairs.open('wb') as out:
        subprocess.check_call(['bgzip', '-c', str(pairs_text)], stdout=out)
    subprocess.check_call(['pairix', '-p', 'pairs', str(pairs)])
    base, mcool = work / 'toy.10.cool', work / 'toy.mcool'
    # Small batches force spill/merge and row-boundary handling in both stages.
    mcool_from_pairs(str(pairs), str(base), str(mcool), '10,20,50', nproc=2,
                     ice_threads=2, window_bp=30, buffer_contacts=7, chunksize=13,
                     io_chunk=5, memory_gb=1, ignore_diags=0, min_nnz=0,
                     min_count=0, mad_max=0, max_iters=1000, tol=1e-8,
                     keep_base_cool=True, temp_dir=str(work / 'scratch'))
    for resolution in (10, 20, 50):
        offsets = {}
        offset = 0
        for name, length in chroms:
            offsets[name] = offset
            offset += (length + resolution - 1) // resolution
        expected = Counter()
        for a, b, x, y, _ in rows:
            i = offsets[a] + (x - 1) // resolution
            j = offsets[b] + (y - 1) // resolution
            expected[(i, j)] += 1
        c = cooler.Cooler(f'{mcool}::resolutions/{resolution}')
        pixels = c.pixels()[:]
        actual = {(int(i), int(j)): int(n) for i, j, n in pixels.itertuples(index=False, name=None)}
        assert actual == dict(expected), f'{resolution}: raw pixels differ from independent counts'
        assert c.info['nbins'] == offset
        assert sum(actual.values()) == len(rows)
        reference, stats = cooler.balance_cooler(c, ignore_diags=0, min_nnz=0, min_count=0,
                                                mad_max=0, tol=1e-8, max_iters=1000, store=False)
        assert stats['converged'], f'{resolution}: reference ICE did not converge'
        weights = c.bins()['weight'][:].to_numpy()
        np.testing.assert_allclose(weights, reference, rtol=1e-8, atol=1e-12, equal_nan=True)
        with h5py.File(mcool, 'r') as f:
            assert bool(f[f'resolutions/{resolution}/bins/weight'].attrs['converged'])
        c.matrix(balance=True, sparse=True).fetch('chr1:0-100')
        print(f'PASS {resolution} bp: raw counts, bins, ICE weights, regional query')
    print('All numerical checks passed. Outputs:', work)


if __name__ == '__main__':
    main()
