# Created on Sun Sep 9 15:04:12 2018

# Author: XiaoTao Wang

import subprocess, os
from runHiC.utilities import sleep

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
            stats[key] = int(value)
        for pre in matchpre:
            if key.startswith(pre):
                stats[key] = int(value)
    
    pipe.communicate()
    
    return stats

def biorep_level(pair_paths, outpre, frag_path):

    # a temporary file to store unfiltered all alignments, for stats and PCR duplicate detection
    out_total = outpre + '.total.pairsam.gz'
    merge_pairs(pair_paths, out_total)
    # cis/trans/dist_freq*/chrom_freq* are all subsets of non duplicates
    # pair_type stats include duplicates
    refkey = ['total', 'total_unmapped', 'total_single_sided_mapped', 'total_mapped', 'total_dups',
              'total_nodups', 'cis', 'trans', 'cis_20kb']
    stats = stats_pairs(out_total, refkey, matchpre=['pair_types','dist_freq'])
    # Final biorep level pairsam
    outpath_1 = outpre + '.pairsam.gz'
    outpath_2 = outpre + '.samfrag.pairsam.gz'

    pipeline = []
    try:
        select_command = ['pairtools', 'select', '\'(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")\'',
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

        select_command = ['pairtools', 'select', '-o', outpath_1, '--output-rest', outpath_2,
                          '\'rfrag_idx1!=rfrag_idx2\'']
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







