Parallel Computation
*********************
*runHiC* provides two layers of parallel computing:

1. You can specify the number of threads / processes when you run ``runHiC mapping`` and
   ``runHiC binning``, by setting the ``-t/--threads`` and the``--nproc`` arguments respectively.
2. If your single machine has limited CPU cores / memory, but you have access to a
   HPC cluster. The most efficient way to process a high-depth sequencing data is to
   submit the same command on multiple nodes. And *runHiC* will parallel the tasks
   for you automatically.

Example
=======
To show how to launch runHiC on our `example data <http://xiaotaowang.github.io/HiC_pipeline/quickstart.html>`_
across multiple nodes of a cluster, let's create a sub-folder and copy the meta data file
into it::

    $ mkdir run-on-cluster
    $ cp datasets.tsv ./run-on-cluster

Then change to the sub-folder. On the 1st node, run command below (you may need to
write a submit script for SLURM/Torque PBS/any other resource manager)::

    $ runHiC pileup -p ../../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 20 -O BAM --include-readid --include-sam --drop-seq --chunkSize 1500000 --logFile runHiC-1.log

Then swith to a 2nd node, and submit the command for the 2nd time::

    $ runHiC pileup -p ../../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 20 -O BAM --include-readid --include-sam --drop-seq --chunkSize 1500000 --logFile runHiC-2.log

Then swith to a 3rd node::

    $ runHiC pileup -p ../../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 20 -O BAM --include-readid --include-sam --drop-seq --chunkSize 1500000 --logFile runHiC-3.log

That's it! You can flexibly allocate the number of threads for each node and the number of
processes to launch to utilize the cluster at the most extent.
