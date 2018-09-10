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

        subprocess.call(' '.join(merge_command), shell=True)

def stats_pairs(inpath, refkey, matchpre):
    
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
    danglingStart = []
    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        pos1 = int(cols[_pairsam_format.COL_P1])
        strand1 = cols[_pairsam_format.COL_S1]
        pos2 = int(cols[_pairsam_format.COL_P2])
        strand2 = cols[_pairsam_format.COL_S2]
        fragstart, fragend = int(cols[-2]), int(cols[-1]) # The index may change in the future
        stats['120_SameFragmentReads'] += 1
        if (strand1=='+') and (strand2=='-'): # dangling reads
            stats['124_DanglingReads'] += 1
            libsize.append(pos2-pos1)
            startsite = min(abs(pos1-fragstart), abs(fragend-pos2))
            danglingStart.append(startsite/(fragend-fragstart))
        elif (strand1=='-') and (strand2=='+'): # self ligation
            stats['122_SelfLigationReads'] += 1
        else:
            stats['126_UnknownMechanism'] += 1
    
    os.remove(samfrag_pairs)
    
    return stats, libsize, danglingStart

def create_frag(chromsizes_file, enzyme):

    # Next


def biorep_level(pair_paths, outpre, frag_path):

    # a temporary file to store unfiltered all alignments, for stats and PCR duplicate detection
    out_total = outpre + '.total.pairsam.gz'
    merge_pairs(pair_paths, out_total)
    # cis/trans/dist_freq*/chrom_freq* are all subsets of non duplicates
    # pair_type stats include duplicates
    refkey = {'total':'000_SequencedReads',
              'total_mapped':'010_DoubleSideMappedReads',
              'total_single_sided_mapped':'020_SingleSideMappedReads',
              'total_unmapped':'030_UnmappedReads',
              'total_dups':'130_DuplicateRemoved',
              'total_nodups':'total_nodups',
              'cis':'410_IntraChromosomalReads',
              'trans':'420_InterChromosomalReads',
              'cis_20kb+':'412_IntraLongRangeReads(>=20Kb)'}

    stats = stats_pairs(out_total, refkey, matchpre=['pair_types','dist_freq'])
    stats['412_IntraShortRangeReads(<20Kb)'] = stats['410_IntraChromosomalReads'] - stats['412_IntraLongRangeReads(>=20Kb)']
    stats['100_NormalPairs'] = stats['010_DoubleSideMappedReads']
    # Final biorep level pairsam
    outpath_1 = outpre + '.pairsam.gz'
    outpath_2 = outpre + '.samfrag.pairsam.gz'

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
        dedup_command = ['pairtools', 'dedup', '--max-mismatch', '2', '--method', 'max']
        pipeline.append(
            subprocess.Popen(dedup_command,
                stdin=pipeline[-1].stdout,
                stdout=subprocess.PIPE,
                bufsize=-1)
        )
        # assign fragment
        restrict_command = ['pairtools', 'restrict', '-f', frag_path]
        pipeline.append(
            subprocess.Popen(restrict_command,
                stdin=pipeline[-1].stdout,
                stdout=subprocess.PIPE,
                bufsize=-1)
        )

        select_command = ['pairtools', 'select', '--output-rest', outpath_1, '-o', outpath_2,
                          '(rfrag_idx1==rfrag_idx2) and (chrom1==chrom2)']
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

    os.remove(out_total)
    os.remove(frag_path)

    substats, libsize, danglingStart = stats_samfrag(outpath_2)

    stats['110_AfterFilteringReads'] = stats['total_nodups'] = substats['120_SameFragmentReads']
    stats['400_TotalContacts'] = stats['110_AfterFilteringReads']
    stats.update(substats)
    del stats['total_nodups']
    # we would never re-count the ligation junction site in original reads
    # we just estimate this ratio
    stats['Ligation Junction Ratio'] = (stats['pair_types/UR'] + stats['pair_types/RU']) / stats['010_DoubleSideMappedReads']
    for key in stats:
        if key.startswith('pair_types'):
            del stats[key]
    
    stats['libsize'] = libsize
    stats['danglingStart'] = danglingStart

    return stats



    

    








