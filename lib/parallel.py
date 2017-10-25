# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:32:51 2017

@author: wxt
"""

import os, pp, subprocess
from collections import Counter

class ppLocal(pp.Server):
    
    def __init__(self, per_worker, maximum_worker):
        
        ncpus = self._detect_ncpus()
        n_worker = ncpus // per_worker
        if not n_worker:
            n_worker = 1
        self.n_worker = min(n_worker, maximum_worker)
        pp.Server.__init__(self, ncpus=self.n_worker)
    
    def _detect_ncpus(self):
        """
        Detects the number of effective CPUs in the system.
        """
        # Linux and Unix
        ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
        if isinstance(ncpus, int) and ncpus > 0:
            return ncpus
        return 1

class ppServer(pp.Server):
    """
    Only support PBS-based clusters.
    """
    def __init__(self, per_worker, maximum_worker, port=60000):
        self.per_worker = per_worker
        self.maximum_worker = maximum_worker
        timeout = self._walltime_to_seconds()
        self._get_nodes()
        servers = self._collect_servers(port)
        self.launch_server(port, timeout)
        pp.Server.__init__(self, ppservers=servers, socket_timeout=timeout)
        
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
    
    def launch_server(self, port, timeout):
        
        template = 'pbsdsh -h {0} ppserver.py -p {1} -w {2} -t {3} -k {4} &'
        self.n_worker = 0
        for node in self.nodes:
            ncpus = self.nodes[node]
            n_worker = ncpus // self.per_worker
            if not n_worker:
                n_worker = 1
            n_worker = min(n_worker, self.maximum_worker)
            self.n_worker += n_worker
            command = template.format(node, port, n_worker, 3600, timeout)
            subprocess.call(command, shell=True)
    
    def _walltime_to_seconds(self):
        
        walltime = os.environ['PBS_WALLTIME']
        if ':' in walltime:
            p = walltime.split(':')
            seconds = 3600*int(p[0]) + 60*int(p[1]) + int(p[2])
        else:
            seconds = int(walltime)
        
        return seconds