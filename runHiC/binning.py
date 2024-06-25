# Created on Thu Sep 13 10:44:19 2018

# Author: XiaoTao Wang

import os, subprocess, logging
from runHiC.utilities import chromsizes_from_pairs

log = logging.getLogger(__name__)

def mcool_from_pairs(pairpath, outcool, outmcool, resolutions, ignore_diags=2, nproc=1, mad_max=5,
                     min_nnz=10, min_count=0, max_split=2):
    
    chromsizes_file, assembly = chromsizes_from_pairs(pairpath)
    
    log.log(21, '{0}: Building contact matrix at a {1} base-pair resolution ...'.format(pairpath, resolutions[0]))
    bin_label = ':'.join([chromsizes_file, str(resolutions[0])])
    bin_command = ['cooler', 'cload', 'pairix', '--assembly', assembly, '--nproc', str(nproc),
                   '--max-split', str(max_split), bin_label, pairpath, outcool]
    subprocess.check_call(' '.join(bin_command), shell=True)

    log.log(21, '{0}: Generating a multi-resolution cooler file by coarsening the {1}bp matrix ...'.format(pairpath, resolutions[0]))
    command = ['cooler', 'zoomify', '-p', str(nproc), '-r {0}'.format(','.join(list(map(str, resolutions)))), '--balance',
               '--balance-args "--mad-max {0} --min-nnz {1} --min-count {2} --ignore-diags {3}"'.format(mad_max, min_nnz, min_count, ignore_diags),
               '-o', outmcool, outcool]
    subprocess.check_call(' '.join(command), shell=True)
    log.log(21, '{0}: Done'.format(pairpath))

    os.remove(outcool)
