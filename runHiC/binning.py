"""Build standard Cooler matrices using the bundled fast runHiC backends."""

import argparse
import json
import logging
import math
import os
import subprocess
import sys

log = logging.getLogger(__name__)
DEFAULT_RESOLUTIONS = '5000,10000,25000,50000,100000,250000,500000,1000000,2500000'


def positive_int(value):
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError('must be a positive integer')
    return result


def nonnegative_int(value):
    result = int(value)
    if result < 0:
        raise argparse.ArgumentTypeError('must be a nonnegative integer')
    return result


def positive_float(value):
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError('must be finite and positive')
    return result


def resolve_resolutions(resolutions, base_resolution=None):
    """Return the loader bin size and the fast_zoomify resolution specification."""
    if isinstance(resolutions, str) and resolutions.upper() in ('B', 'N'):
        if base_resolution is None or base_resolution <= 0:
            raise ValueError('--base-resolution is required with --resolutions B or N')
        return base_resolution, resolutions.upper()
    try:
        values = sorted(set(int(r) for r in (
            resolutions.split(',') if isinstance(resolutions, str) else resolutions)))
    except (TypeError, ValueError):
        raise ValueError('--resolutions must be B, N, or comma-separated positive integers')
    if not values or values[0] <= 0:
        raise ValueError('Resolutions must be positive integers')
    base = values[0] if base_resolution is None else base_resolution
    if base <= 0 or any(r < base or r % base for r in values):
        raise ValueError('Every resolution must be an integer multiple of the base resolution')
    return base, ','.join(map(str, sorted(set([base] + values))))


def add_binning_arguments(parser):
    """Use the same matrix options for binning and the full pileup pipeline."""
    parser.add_argument('--resolutions', default=DEFAULT_RESOLUTIONS,
                        help='Comma-separated bp sizes, or B/N with --base-resolution. Base always included. Defaults retain the original 5 kb pyramid; specify 10 for 10 bp.')
    parser.add_argument('--base-resolution', type=positive_int, help='Loader bin size in bp; otherwise the smallest explicit resolution.')
    parser.add_argument('--nproc', type=positive_int, default=8, help='Pairix loader workers (also filtering workers in pileup).')
    parser.add_argument('--ice-threads', type=positive_int, help='ICE threads; defaults to --nproc. Coarsening is single-process.')
    parser.add_argument('--max-split', type=positive_int, help='Deprecated and ignored; tune --window-bp and --buffer-contacts instead.')
    parser.add_argument('--window-bp', type=positive_int, default=1000000, help='Pairix query interval width in bp.')
    parser.add_argument('--buffer-contacts', type=positive_int, default=2000000, help='Contacts per worker buffer before aggregation/spilling.')
    parser.add_argument('--chunksize', type=positive_int, default=2000000, help='Pixels per coarsening/ICE batch; distinct from mapping --chunkSize.')
    parser.add_argument('--io-chunk', type=positive_int, default=250000, help='HDF5/index/spill-merge batch length for both backends.')
    parser.add_argument('--compression', type=int, choices=range(10), default=1, help='Gzip level for newly written matrices and weights.')
    parser.add_argument('--memory-gb', type=positive_float, default=128, help='ICE planning budget in GiB, not a hard cap. Distinct from pairtools --memory.')
    parser.add_argument('--temp-dir', help='Matrix scratch; defaults to system temp. Distinct from mapping/filtering --tmpdir.')
    parser.add_argument('-D', '--ignore-diags', type=nonnegative_int, default=2, help='Diagonals excluded from ICE (0 includes all).')
    parser.add_argument('--ignore-dist', type=nonnegative_int, default=0, help='ICE exclusion distance in bp, combined with --ignore-diags.')
    parser.add_argument('--min-nnz', type=nonnegative_int, default=10, help='Minimum marginal nonzeros before ICE.')
    parser.add_argument('--min-count', type=nonnegative_int, default=0, help='Minimum marginal counts before ICE.')
    parser.add_argument('--mad-max', type=nonnegative_int, default=5, help='MAD filter threshold before ICE.')
    parser.add_argument('--blacklist', help='Headerless BED of regions to exclude from ICE.')
    parser.add_argument('--tol', type=positive_float, default=1e-5, help='ICE marginal-variance convergence tolerance.')
    parser.add_argument('--max-iters', type=positive_int, default=200, help='Maximum ICE iterations per resolution.')
    parser.add_argument('--convergence-policy', choices=('store_nan', 'error', 'store_final'), default='store_nan',
                        help='On nonconvergence: write NaNs, retain finite final weights at the iteration limit, or abort without writing weights. All policies report stage failure.')
    parser.add_argument('--resume', action='store_true', help='Resume a fast mcool with its retained base cool and provenance. New samples start normally.')
    parser.add_argument('--keep-base-cool', action='store_true', help='Keep the intermediate base cool and provenance even after success.')


def _fingerprint(path):
    stat = os.stat(path)
    return dict(path=os.path.realpath(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def _run_backend(command):
    # Argument lists preserve spaces/non-ASCII paths without invoking a shell.
    log.log(21, 'Running matrix backend: %s', subprocess.list2cmdline(command))
    subprocess.check_call(command)


def mcool_from_pairs(pairpath, outcool, outmcool, resolutions, ignore_diags=2, nproc=1,
                     mad_max=5, min_nnz=10, min_count=0, max_split=None,
                     base_resolution=None, ice_threads=None, window_bp=1000000,
                     buffer_contacts=2000000, chunksize=2000000, io_chunk=250000,
                     compression=1, memory_gb=128, temp_dir=None, ignore_dist=0,
                     blacklist=None, tol=1e-5, max_iters=200,
                     convergence_policy='store_nan', resume=False, keep_base_cool=False):
    if sys.version_info < (3, 10):
        raise RuntimeError('The fast matrix backends require Python 3.10 or newer')
    base, resolution_spec = resolve_resolutions(resolutions, base_resolution)
    if max_split is not None:
        log.warning('--max-split is ignored by the fast loader; use --window-bp and --buffer-contacts')
    pairpath, outcool, outmcool = map(os.path.abspath, (pairpath, outcool, outmcool))
    if len(set(map(os.path.normcase, (pairpath, outcool, outmcool)))) != 3:
        raise ValueError('Pairs, base .cool and output .mcool must be different paths')
    index = pairpath + '.px2'
    if not os.path.isfile(index):
        raise FileNotFoundError('Pairix index is required: ' + index)
    provenance = dict(pairs=_fingerprint(pairpath), index=_fingerprint(index),
                      base_resolution=base, coordinates='one-based')
    provenance_path = outcool + '.runhic.json'
    if os.path.exists(outmcool) or (resume and os.path.exists(outcool)):
        if not resume:
            raise FileExistsError('Output exists; use --resume for a matching fast run or choose a new output: ' + outmcool)
        if not os.path.isfile(outcool) or not os.path.isfile(provenance_path):
            raise FileNotFoundError('Resume requires the retained base .cool and its .runhic.json provenance')
        with open(provenance_path) as stream:
            saved = json.load(stream)
        if saved != dict(source=provenance, base=_fingerprint(outcool)):
            raise ValueError('Pairs/index/base changed since loading; choose a new output instead of resuming')
        log.log(21, 'Resuming %s with retained base %s', outmcool, outcool)
    else:
        if os.path.exists(outcool) or os.path.exists(provenance_path):
            raise FileExistsError('Base output already exists without an mcool to resume; move it aside or use fast_zoomify directly: ' + outcool)
        from runHiC.utilities import chromsizes_from_pairs
        chromsizes_file, assembly = chromsizes_from_pairs(pairpath)
        log.log(21, '%s: Building the %d bp matrix with fast_pairix_to_cool', pairpath, base)
        command = [sys.executable, os.path.join(os.path.dirname(__file__), 'fast_pairix_to_cool.py'),
                   '--assembly', assembly, '--nproc', str(nproc),
                   '--window-bp', str(window_bp), '--buffer-contacts', str(buffer_contacts),
                   '--io-chunk', str(io_chunk), '--compression', str(compression)]
        if temp_dir:
            command.extend(['--temp-dir', os.path.abspath(os.path.expanduser(temp_dir))])
        command.extend([chromsizes_file + ':' + str(base), pairpath, outcool])
        _run_backend(command)
        if provenance != dict(pairs=_fingerprint(pairpath), index=_fingerprint(index),
                              base_resolution=base, coordinates='one-based'):
            raise RuntimeError('Pairs or index changed while loading; base matrix cannot be reused safely')
        with open(provenance_path, 'x') as stream:
            json.dump(dict(source=provenance, base=_fingerprint(outcool)), stream, sort_keys=True)

    log.log(21, '%s: Coarsening and balancing every requested resolution with fast_zoomify', pairpath)
    command = [sys.executable, os.path.join(os.path.dirname(__file__), 'fast_zoomify.py'), outcool,
               '--out', outmcool, '--resolutions', resolution_spec,
               '--threads', str(nproc if ice_threads is None else ice_threads),
               '--chunksize', str(chunksize), '--io-chunk', str(io_chunk),
               '--compression', str(compression), '--memory-gb', str(memory_gb),
               '--ignore-diags', str(ignore_diags), '--ignore-dist', str(ignore_dist),
               '--min-nnz', str(min_nnz), '--min-count', str(min_count), '--mad-max', str(mad_max),
               '--tol', str(tol), '--max-iters', str(max_iters), '--convergence-policy', convergence_policy]
    if temp_dir:
        command.extend(['--temp-dir', os.path.abspath(os.path.expanduser(temp_dir))])
    if blacklist:
        command.extend(['--blacklist', os.path.abspath(os.path.expanduser(blacklist))])
    if resume:
        command.append('--resume')
    try:
        _run_backend(command)
    except subprocess.CalledProcessError as exc:
        if exc.returncode == 2:
            log.error('ICE did not converge at all levels. Raw mcool/base matrices are retained; no completion marker will be written. Inspect weights and retry with --resume.')
        raise
    if not keep_base_cool:
        os.remove(outcool)
        os.remove(provenance_path)
    log.log(21, '%s: Done', pairpath)
