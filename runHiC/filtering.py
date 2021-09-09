# Created on Sun Sep 9 15:04:12 2018

# Author: XiaoTao Wang

import subprocess, os, pickle
from runHiC.utilities import sleep
import numpy as np
from copy import deepcopy
import time
from collections import defaultdict

def merge_pairs(pair_paths, outpath, tmpdir, nproc_in, nproc_out, memory):

    if len(pair_paths)==1:
        # Just make a soft link to original .pairsam
        os.symlink(pair_paths[0], outpath)
    else:
        # runHiC doesn't provide interface for changing detailed parameters of
        # pairtools merge for simplicity
        merge_command = ['pairtools', 'merge', '-o', outpath, '--nproc', str(nproc_out), '--memory', memory,
                         '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out),
                         '--max-nmerge', '8', '--tmpdir', tmpdir] + pair_paths
        subprocess.check_call(' '.join(merge_command), shell=True)


def dedup(out_total, outpath, stats, nproc_in, nproc_out):

    pipeline = []
    try:
        dedup_command = ['pairtools', 'dedup', '--max-mismatch', '1', '--method', 'max',
                         '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out),
                         '-o', outpath, out_total]
        pipeline.append(
            subprocess.Popen(dedup_command,
                stdout=None,
                bufsize=-1)
        )
        pipeline[-1].wait()
    finally:
        sleep()
        for process in pipeline:
            if process.poll() is None:
                process.terminate()
    
    os.remove(out_total)

    refkey = {'cis':'410_IntraChromosomalReads',
              'trans':'420_InterChromosomalReads',
              'cis_20kb+':'412_IntraLongRangeReads(>=20Kb)',
              'total_nodups':'total_nodups'}
    
    substats = stats_pairs(outpath, refkey, matchpre=['dist_freq'], nproc_in=nproc_in, nproc_out=nproc_out)
    stats['130_DuplicateRemoved'] = stats['110_AfterFilteringReads'] - substats['total_nodups']
    stats['110_AfterFilteringReads'] = substats['total_nodups']
    stats['400_TotalContacts'] = stats['110_AfterFilteringReads']
    stats.update(substats)
    stats['412_IntraShortRangeReads(<20Kb)'] = stats['410_IntraChromosomalReads'] - stats['412_IntraLongRangeReads(>=20Kb)']
    del stats['total_nodups']


def collect_stats(pair_paths):

    stats_pool = {}
    for i, p in enumerate(pair_paths):
        stats_path = p.replace('.pairsam.gz', '.pstats.1')
        with open(stats_path, 'rb') as source:
            stats_pool[str(i)] = pickle.load(source)['pseudo']
    
    merge_stats(stats_pool, list(stats_pool.keys()), 'pseudo')
    stats_pool = {'pseudo': stats_pool['pseudo']}

    return stats_pool


def stats_pairs(inpath, refkey, matchpre=[], nproc_in=3, nproc_out=8):
    
    stat_command = ['pairtools', 'stats', '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out), inpath]
    pipe = subprocess.Popen(stat_command, stdout=subprocess.PIPE)
    inStream = pipe.stdout
    stats = defaultdict(int)
    for line in inStream:
        parse = line.decode().rstrip().split()
        key, value = parse[0], parse[1]
        if key in refkey:
            stats[refkey[key]] = int(value)
        for pre in matchpre:
            if key.startswith(pre):
                stats[key] = int(value)
    
    pipe.communicate()
    
    return stats

def stats_samfrag(samfrag_pairs, sample_size=100000):

    from pairtools import _fileio, _pairsam_format, _headerops

    instream = _fileio.auto_open(samfrag_pairs, mode='r')
    _, body_stream = _headerops.get_header(instream)

    stats = defaultdict(int)
    libsize = []
    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        pos1 = int(cols[_pairsam_format.COL_P1])
        strand1 = cols[_pairsam_format.COL_S1]
        pos2 = int(cols[_pairsam_format.COL_P2])
        strand2 = cols[_pairsam_format.COL_S2]
        #fragstart, fragend = int(cols[-2]), int(cols[-1]) # The index may change in the future
        stats['120_SameFragmentReads'] += 1
        if (strand1=='+') and (strand2=='-'): # dangling reads
            stats['124_DanglingReads'] += 1
            libsize.append(pos2-pos1)
        elif (strand1=='-') and (strand2=='+'): # self ligation
            stats['122_SelfLigationReads'] += 1
        else:
            stats['126_UnknownMechanism'] += 1
    
    instream.close()
    
    os.remove(samfrag_pairs)

    libsize = np.r_[libsize]
    np.random.shuffle(libsize)
    libsize = libsize[:sample_size]
    
    return stats, libsize

def create_frag(genomepath, chromsizes_file, enzyme, tmpdir):

    _, fastaName = os.path.split(genomepath)
    genomeName = fastaName.split('.fa')[0]
    prefix = '.'.join([genomeName, 'frags', enzyme])
    outbed = os.path.join(tmpdir, '.'.join([prefix, 'bed']))
    lockfil = os.path.join(tmpdir, '.'.join([prefix, 'lock']))
    if os.path.exists(outbed) or os.path.exists(lockfil):
        while os.path.exists(lockfil):
            time.sleep(0.5)
        return outbed
    else:
        lock = open(lockfil, 'wb') # aquire the file lock
        lock.close()
        digest_command = ['runHiC-digest', '-O', outbed, '-C', chromsizes_file, '--fasta-path', genomepath, '--enzyme', enzyme]
        subprocess.check_call(' '.join(digest_command), shell=True)
        os.remove(lockfil)

        return outbed

def biorep_level(pair_paths, outpre, tmpdir, nproc_in, nproc_out, memory):

     # Final biorep level pairsam
    outpath = outpre + '.pairsam.gz'
    merge_pairs(pair_paths, outpath, tmpdir, nproc_in, nproc_out, memory)
    stats = collect_stats(pair_paths)['pseudo']
    
    return stats, outpath

def merge_stats(stats_pool, keys, outkey, sample_size=100000):

    stats_pool[outkey] = deepcopy(stats_pool[keys[0]])
    for i in stats_pool[outkey].keys():
        for k in keys[1:]:
            if not i in stats_pool[k]:
                continue
            if not i in ['libsize']:
                stats_pool[outkey][i] += stats_pool[k][i]
            else:
                stats_pool[outkey][i] = np.r_[stats_pool[outkey][i], stats_pool[k][i]]
    
    if 'libsize' in stats_pool[outkey]:
        np.random.shuffle(stats_pool[outkey]['libsize'])
        stats_pool[outkey]['libsize'] = stats_pool[outkey]['libsize'][:sample_size] # limit sample size
    
def enzyme_level(pair_paths, outpre, keys, outkey, stats_pool, tmpdir, nproc_in, nproc_out, memory):

    ## pair_paths --> outpre
    ## keys --> outkey
    outall = outpre + '.pairsam.gz'
    merge_pairs(pair_paths, outall, tmpdir, nproc_in, nproc_out, memory)
    merge_stats(stats_pool, keys, outkey)

    return stats_pool, outall

def split_pairsam(pairsam_path):

    from pairtools import _fileio, _headerops
    
    # check for SAM information with the header of .pairsam.gz
    instream = _fileio.auto_open(pairsam_path, mode='r')
    header, _ = _headerops.get_header(instream)
    SAM = False
    for r in header:
        if not r.startswith('#columns:'):
            continue
        columns = r.split(':')[1].strip().split()
        if ('sam1' in columns) and ('sam2' in columns):
            SAM = True
    
    instream.close()
    pairpath = pairsam_path.replace('.pairsam.gz', '.pairs.gz')
    if SAM:
        bampath = pairsam_path.replace('.pairsam.gz', '.bam')
        split_command = ['pairtools', 'split', '--output-pairs', pairpath,
                         '--output-sam', bampath, pairsam_path]
        subprocess.check_call(' '.join(split_command), shell=True)
    else:
        mv_command = ['cp', pairsam_path, pairpath]
        subprocess.check_call(' '.join(mv_command), shell=True)
    
    # generate pairix index
    pairix_command = ['pairix', pairpath]
    subprocess.check_call(' '.join(pairix_command), shell=True)

    return pairpath



