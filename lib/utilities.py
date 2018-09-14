# Created on Tue Dec 23 21:15:19 2014

# Author: XiaoTao Wang

import logging, os, subprocess, time, gc
import numpy as np
from mirnylib.h5dict import h5dict

log = logging.getLogger(__name__)

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

def chromsizes_from_fasta(genomeFolder, genomeName):

    from Bio import SeqIO

    genomepath = os.path.join(genomeFolder, '.'.join([genomeName, 'fa']))
    chromsizes = {}
    fasta_genome = SeqIO.parse(genomepath, 'fasta')
    for fasta in fasta_genome:
        chromsizes[fasta.id] = len(fasta.seq)
    
    outfile = os.path.join(genomeFolder, '.'.join([genomeName, 'chrom', 'sizes']))
    with open(outfile, 'w') as out:
        for c in sorted(chromsizes):
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