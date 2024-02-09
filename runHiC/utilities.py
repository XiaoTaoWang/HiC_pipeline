# Created on Tue Dec 23 21:15:19 2014

# Author: XiaoTao Wang

import os, time, gc, subprocess, tempfile
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

    import pairtools

    if pairtools.__version__.startswith('0'):
        from pairtools import _fileio, _headerops
    else:
        from pairtools.lib import fileio as _fileio
        from pairtools.lib import headerops as _headerops

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

def chrname_sort_flip_order(genomepath, chrom_fil, tmpdir):

    pre_defined = [line.rstrip().split()[0] for line in open(chrom_fil, 'r')]
    ref_chroms = set()
    with open(genomepath, 'r') as source:
        for line in source:
            if line.startswith('>'):
                chrom_name = line.rstrip().lstrip('>')
                ref_chroms.add(chrom_name)
    
    tl = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
    kw = {'suffix':tl, 'dir':tmpdir}
    fd, flip_order_fil = tempfile.mkstemp(**kw)
    os.close(fd)
    flip_order = []
    with open(flip_order_fil, 'w') as out:
        for c in pre_defined:
            if c in ref_chroms:
                out.write(c+'\n')
                flip_order.append(c)
    
    fd, sort_order_fil = tempfile.mkstemp(**kw)
    os.close(fd)
    with open(sort_order_fil, 'w') as out:
        for c in sorted(flip_order):
            out.write(c+'\n')
    
    return flip_order_fil, sort_order_fil

def reorder_chromosomes_in_fasta(genomepath, chrom_fil, tmpdir):

    chrom_order = [line.rstrip().split()[0] for line in open(chrom_fil, 'r')]
    # extract reference sequences by chromosome
    outstream = None
    ref_chroms = set()
    with open(genomepath, 'r') as source:
        for line in source:
            if line.startswith('>'):
                if not outstream is None:
                    outstream.close()
                chrom_name = line.rstrip().lstrip('>')
                outstream = open(os.path.join(tmpdir, '{0}.fa'.format(chrom_name)), 'w')
                outstream.write(line)
                ref_chroms.add(chrom_name)
            else:
                if not outstream is None:
                    outstream.write(line)
    
    if not outstream is None:
        outstream.close()
    
    # out chromosome order
    chroms = []
    for c in chrom_order:
        if c in ref_chroms:
            chroms.append(c)
    # output
    command = ['cat'] + [os.path.join(tmpdir, '{0}.fa'.format(c)) for c in chroms] + ['>', genomepath]
    subprocess.check_call(' '.join(command), shell=True)
    for c in ref_chroms:
        os.remove(os.path.join(tmpdir, '{0}.fa'.format(c)))


def extract_chrom_sizes(fil):

    chromsizes = []
    with open(fil, 'r') as source:
        for line in source:
            tmp = line.rstrip().split()
            chromsizes.append((tmp[0], int(tmp[1])))
    
    return chromsizes

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
            if '-' in enzyme:
                cocktail = biorst.RestrictionBatch(enzyme.split('-'))
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