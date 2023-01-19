Parallel Computation
*********************
*runHiC* provides two layers of parallel computing:

1. You can specify the number of threads / processes when you run ``runHiC mapping``,
   ``runHiC filtering``, ``runHiC binning`` and ``runHiC pileup`` subcommands, by setting
   the ``-t/--threads`` or ``--nproc`` arguments.
2. If you have access to a HPC cluster. The most efficient way to process a
   deeply-sequenced library is to submit the same command across multiple nodes.
   And *runHiC* will parallel the tasks for you automatically.

Example
=======
To show how to launch runHiC on our `example data <http://xiaotaowang.github.io/HiC_pipeline/quickstart.html>`_
across multiple nodes of a cluster, let's create a sub-folder within *workspace* and copy the meta data file
into it::

    $ mkdir run-on-cluster
    $ cp datasets.tsv ./run-on-cluster
    $ cd run-on-cluster

Then change to the *run-on-cluster* sub-folder, and submit the command below (you may need
to write a job submission script based on the resource management system installed on your
cluster)::

    $ runHiC pileup -p ../../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 20 --include-readid --drop-seq --chunkSize 1500000 --logFile runHiC-1.log

Then swith to a different node, and submit the same command for the 2nd time::

    $ runHiC pileup -p ../../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 20 --include-readid --drop-seq --chunkSize 1500000 --logFile runHiC-2.log

Then swith to a 3rd node::

    $ runHiC pileup -p ../../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 20 --include-readid --drop-seq --chunkSize 1500000 --logFile runHiC-3.log

That's it! You can flexibly allocate the number of threads for each node and the number of
processes to launch to utilize the cluster at the most extent.
