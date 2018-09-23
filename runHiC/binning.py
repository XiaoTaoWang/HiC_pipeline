# Created on Thu Sep 13 10:44:19 2018

# Author: XiaoTao Wang

import os, subprocess, logging
from runHiC.utilities import chromsizes_from_pairs

log = logging.getLogger(__name__)

def binning_from_pairs(pairpath, res, outpath, ignore_diags=1, nproc=1, mad_max=5, min_nnz=10,
                       min_count=0):

    chromsizes_file, assembly = chromsizes_from_pairs(pairpath)
    
    log.log(21, 'Building contact matrix ...')
    bin_label = ':'.join([chromsizes_file, str(res)])
    bin_command = ['cooler', 'cload', 'pairix', '--assembly', assembly, '--nproc', str(nproc),
                   bin_label, pairpath, outpath]
    subprocess.check_call(' '.join(bin_command), shell=True)
    log.log(21, 'Done')

    log.log(21, 'Perform ICE ...')
    ice_command = ['cooler', 'balance', '--nproc', str(nproc), '--ignore-diags', str(ignore_diags),
                   '--mad-max', str(mad_max), '--min-nnz', str(min_nnz), '--min-count', str(min_count),
                   outpath]
    subprocess.check_call(' '.join(ice_command), shell=True)
    log.log(21, 'Done')



