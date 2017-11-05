# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:32:51 2017

@author: wxt
"""

import os, time, pp, random, subprocess, logging
import multiprocessing as mp
from collections import Counter

log = logging.getLogger(__name__)

def mpPool(per_worker, maximum_worker):
    
    local = os.environ["HOSTNAME"]
    ncpus = mp.cpu_count()
    n_worker = ncpus // per_worker
    if not n_worker:
        n_worker = 1
    n_worker = min(n_worker, maximum_worker)
    log.log(21, 'Launch %d processes on local compute node: %s',
            n_worker, local)
    
    return mp.Pool(n_worker), n_worker
        

class ppLocal(pp.Server):
    
    def __init__(self, per_worker, maximum_worker):
        
        local = os.environ["HOSTNAME"]
        ncpus = mp.cpu_count()
        n_worker = ncpus // per_worker
        if not n_worker:
            n_worker = 1
        self.n_worker = min(n_worker, maximum_worker)
        log.log(21, 'Launch %d processes on local compute node: %s',
                self.n_worker, local)
        pp.Server.__init__(self, ncpus=self.n_worker)
        

class ppServer(pp.Server):
    """
    Only support PBS-based clusters.
    """
    def __init__(self, per_worker, maximum_worker, port=60000):
        self.per_worker = per_worker
        self.maximum_worker = maximum_worker
        secret = self._genSecret()
        timeout = self._walltime_to_seconds()
        self._get_nodes()
        #local_worker = self._local_node()
        servers = self._collect_servers(port)
        self.launch_server(port, secret, timeout)
        pp.Server.__init__(self, ncpus=0, ppservers=servers, secret=secret,
                           socket_timeout=timeout)
    
    def _cal_worker(self, ncpus):
        
        n_worker = ncpus // self.per_worker
        if not n_worker:
            n_worker = 1
        n_worker = min(n_worker, self.maximum_worker)
        
        return n_worker
    
    def _local_node(self):
        
        local = os.environ["HOSTNAME"]
        local_ncpus = self.nodes[local]
        self.nodes.pop(local)
        n_worker = self._cal_worker(local_ncpus)
        
        log.log(21, 'Launch %d processes on local compute node: %s',
                n_worker, local)
        
        return n_worker
        
    def _get_nodes(self):
        pbs_nodefile = os.environ["PBS_NODEFILE"]
        procs = [line.rstrip() for line in open(pbs_nodefile,'r')]
        self.nodes = Counter(procs)
    
    def _collect_servers(self, port):
        
        servers = ()
        for node in self.nodes:
            tmp = '{0}:{1}'.format(node, port)
            servers += (tmp,)
        return servers
    
    def launch_server(self, port, secret, timeout):
        
        template = 'pbsdsh -h {0} ppserver.py -p {1} -w {2} -s {3} -t {4} -k {5} &'
        for node in self.nodes:
            ncpus = self.nodes[node]
            n_worker = self._cal_worker(ncpus)
            log.log(21, 'Launch %d processed on remote compute node: %s',
                    n_worker, node)
            command = template.format(node, port, n_worker, secret, 3600, timeout)
            subprocess.call(command, shell=True)
            time.sleep(5)
    
    def _genSecret(self):
        
        tl = list(time.strftime('%Y%m%d%H%M%S',time.localtime(time.time())))
        random.shuffle(tl)
        secret = ''.join(tl)

        return secret
    
    def _walltime_to_seconds(self):
        
        walltime = os.environ['PBS_WALLTIME']
        if ':' in walltime:
            p = walltime.split(':')
            seconds = 3600*int(p[0]) + 60*int(p[1]) + int(p[2])
        else:
            seconds = int(walltime)
        
        return seconds