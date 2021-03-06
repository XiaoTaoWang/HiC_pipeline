# Created on Tue Dec 23 21:15:19 2014

# Author: XiaoTao Wang

import os, time, gc
import numpy as np
import pandas as pd
from cooler.util import load_fasta, read_chromsizes

def sleep():
    for _ in range(3):
        time.sleep(0.1)
    gc.collect()
    for _ in range(3):
        time.sleep(0.1)

def cleanFile(filename):
    if os.path.exists(filename):
        os.remove(filename)
    
def cleanDirectory(dirName):
    for i in os.listdir(dirName):
        os.remove(os.path.join(dirName, i))

def find_digit_parts(chrname):

    collect = []
    for s in chrname[::-1]:
        if s.isdigit():
            collect.append(s)
        else:
            break
    
    if len(collect):
        digit_parts = int(''.join(collect[::-1]))
        return digit_parts
    else:
        return

def sort_chromlabels(chrnames):

    num_table = []
    nonnum_names = []
    for n in chrnames:
        tmp = find_digit_parts(n)
        if tmp is None:
            nonnum_names.append(n)
        else:
            num_table.append((tmp, n))

    num_table.sort()
    sorted_names = [s[1] for s in num_table]

    for s in ['M', 'Y', 'X']:
        for idx, n in enumerate(nonnum_names):
            if n.endswith(s):
                nonnum_names.pop(idx)
                nonnum_names.insert(0, n)
                break
    sorted_names = sorted_names + nonnum_names

    return sorted_names


def chromsizes_from_fasta(genomeFolder, genomeName):

    from Bio import SeqIO

    genomepath = os.path.join(genomeFolder, '.'.join([genomeName, 'fa']))
    chromsizes = {}
    fasta_genome = SeqIO.parse(genomepath, 'fasta')
    for fasta in fasta_genome:
        chromsizes[fasta.id] = len(fasta.seq)
    
    outfile = os.path.join(genomeFolder, '.'.join([genomeName, 'chrom', 'sizes']))
    with open(outfile, 'w') as out:
        for c in sort_chromlabels(chromsizes):
            line = [c, str(chromsizes[c])]
            out.write('\t'.join(line)+'\n')
    
    return outfile

def chromsizes_from_pairs(pairpath):

    from pairtools import _fileio, _headerops

    instream = _fileio.auto_open(pairpath, mode='r')
    header, _ = _headerops.get_header(instream)
    genomeName = 'Unknown'
    chromsizes = []
    for r in header:
        if r.startswith('#genome_assembly:'):
            genomeName = r.split(':')[1].strip()
        if r.startswith('#chromsize:'):
            pair = r.split(':')[1].strip().split()
            chromsizes.append(pair)
    
    folder = os.path.split(pairpath)[0]
    outpath = os.path.join(folder, '.'+genomeName+'.chrom.sizes') # invisible to users
    with open(outpath, 'w') as out:
        for line in chromsizes: # order unchanged
            out.write('\t'.join(line)+'\n')
    
    instream.close()
    
    return outpath, genomeName

def digest(fasta_records, enzyme):
    """
    Divide a genome into restriction fragments. Support Arima-HiC enzyme cocktail
    which digest chromatin at ^GATC and G^ANTC.

    Parameters
    ----------
    fasta_records : OrderedDict
        Dictionary of chromosome names to sequence records.

    enzyme: str
        Name of restriction enzyme.

    Returns
    -------
    Dataframe with columns: 'chrom', 'start', 'end'.

    """
    import Bio.Restriction as biorst
    import Bio.Seq as bioseq
    # http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId447698
    chroms = fasta_records.keys()
    try:
        if enzyme.lower()=='arima':
            cocktail = biorst.RestrictionBatch(['MboI', 'HinfI'])
            cut_finder = cocktail.search
        else:
            cut_finder = getattr(biorst, enzyme).search
    except AttributeError:
        raise ValueError('Unknown enzyme name: {}'.format(enzyme))

    def _each(chrom):
        seq = bioseq.Seq(str(fasta_records[chrom]))
        tmp = cut_finder(seq)
        if type(tmp)==list:
            cut_sites = tmp
        elif type(tmp)==dict:
            cut_sites = []
            for e in tmp:
                cut_sites.extend(tmp[e])
            cut_sites.sort()
        cuts = np.r_[0, np.array(cut_sites) + 1, len(seq)].astype(int)
        n_frags = len(cuts) - 1

        frags = pd.DataFrame({
            'chrom': [chrom] * n_frags,
            'start': cuts[:-1],
            'end': cuts[1:]},
            columns=['chrom', 'start', 'end'])
        return frags

    return pd.concat(map(_each, chroms), axis=0, ignore_index=True)