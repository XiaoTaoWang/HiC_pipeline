# Created on Sun Sep 9 15:04:12 2018

# Author: XiaoTao Wang

import subprocess, os
from runHiC.utilities import sleep
import numpy as np

def merge_pairs(pair_paths, outpath):

    if len(pair_paths)==1:
        # Just make a soft link to original .pairsam
        os.symlink(pair_paths[0], outpath)
    else:
        # runHiC doesn't provide interface for changing detailed parameters of
        # pairtools merge for simplicity
        merge_command = ['pairtools', 'merge', '-o', outpath, '--nproc', '8', '--memory', '2G',
                         '--max-nmerge', '8'] + pair_paths

        subprocess.check_call(' '.join(merge_command), shell=True)

def stats_pairs(inpath, refkey, matchpre=[]):
    
    stat_command = ['pairtools', 'stats', inpath]
    pipe = subprocess.Popen(stat_command, stdout=subprocess.PIPE)
    inStream = pipe.stdout
    stats = {}
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

def stats_samfrag(samfrag_pairs):

    from pairtools import _fileio, _pairsam_format, _headerops
    from collections import defaultdict

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
    
    return stats, libsize

def create_frag(genomepath, chromsizes_file, enzyme):

    Folder, fastaName = os.path.split(genomepath)
    genomeName = fastaName.split('.')[0]
    outbed = os.path.join(Folder, '.'.join([genomeName, 'frags', enzyme, 'bed']))
    if os.path.exists(outbed):
        return outbed
    
    digest_command = ['cooler', 'digest', '-o', outbed, chromsizes_file, genomepath, enzyme]
    subprocess.check_call(' '.join(digest_command), shell=True)

    return outbed

def biorep_level(pair_paths, outpre, frag_path):

    # a temporary file to store unfiltered all alignments
    out_total = outpre + '.total.pairsam.gz'
    merge_pairs(pair_paths, out_total)
    # pair_type stats include duplicates
    refkey = {'total':'000_SequencedReads',
              'total_mapped':'010_DoubleSideMappedReads',
              'total_single_sided_mapped':'020_SingleSideMappedReads',
              'total_unmapped':'030_UnmappedReads'
              }

    stats = stats_pairs(out_total, refkey)
    stats['100_NormalPairs'] = stats['010_DoubleSideMappedReads']
    # Final biorep level pairsam
    outpath = outpre + '.pairsam.gz' # select.dedup.filter
    outpath_1 = outpre + '.select.dedup.pairsam.gz'
    outpath_2 = outpre + '.select.dedup.samefrag.pairsam.gz'

    pipeline = []
    try:
        select_command = ['pairtools', 'select', '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")',
                          out_total]
        pipeline.append(
            subprocess.Popen(select_command,
                stdout=subprocess.PIPE,
                bufsize=-1
            )
        )
        dedup_command = ['pairtools', 'dedup', '--max-mismatch', '1', '--method', 'max', '-o', outpath_1]
        pipeline.append(
            subprocess.Popen(dedup_command,
                stdin=pipeline[-1].stdout,
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
    
    refkey = {'total_nodups':'total_nodups'}
    dupstats = stats_pairs(outpath_1, refkey)
    stats['130_DuplicateRemoved'] = stats['100_NormalPairs'] - dupstats['total_nodups']

    try:
        # assign fragment
        restrict_command = ['pairtools', 'restrict', '-f', frag_path, outpath_1]
        pipeline.append(
            subprocess.Popen(restrict_command,
                stdout=subprocess.PIPE,
                bufsize=-1)
        )

        ####### COLS[-6]==COLS[-3], the index may change to follow pairtools
        select_command = ['pairtools', 'select', '--output-rest', outpath, '-o', outpath_2,
                          '(COLS[-6]==COLS[-3]) and (chrom1==chrom2)']
        pipeline.append(
            subprocess.Popen(select_command,
                stdin=pipeline[-1].stdout,
                stdout=None,
                bufsize=-1)
        )
        pipeline[-1].wait()
    finally:
        sleep()
        for process in pipeline:
            if process.poll() is None:
                process.terminate()

    os.remove(outpath_1)

    substats, libsize = stats_samfrag(outpath_2)

    stats['110_AfterFilteringReads'] = dupstats['total_nodups'] - substats['120_SameFragmentReads']
    stats['400_TotalContacts'] = stats['110_AfterFilteringReads']
    stats.update(substats)

    refkey = {'cis':'410_IntraChromosomalReads',
              'trans':'420_InterChromosomalReads',
              'cis_20kb+':'412_IntraLongRangeReads(>=20Kb)'
              }
    
    substats = stats_pairs(outpath, refkey, matchpre=['dist_freq'])
    stats.update(substats)
    stats['412_IntraShortRangeReads(<20Kb)'] = stats['410_IntraChromosomalReads'] - stats['412_IntraLongRangeReads(>=20Kb)']

    # we would never re-count the ligation junction site in original reads
    
    stats['libsize'] = libsize

    return stats, outpath
    

def enzyme_level(pair_paths, outpre, keys, outkey, stats_pool):

    ## pair_paths --> outpre
    ## keys --> outkey
    outall = outpre + '.pairsam.gz'
    merge_pairs(pair_paths, outall)
    stats_pool[outkey] = stats_pool[keys[0]]
    for i in stats_pool[outkey].keys():
        for k in keys[1:]:
            if not i in ['libsize']:
                stats_pool[outkey][i] += stats_pool[k][i]
            else:
                stats_pool[outkey][i] = np.r_[stats_pool[outkey][i], stats_pool[k][i]]

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



