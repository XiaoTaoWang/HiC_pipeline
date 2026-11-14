Parallel Computation
********************

*runHiC* provides two layers of parallel computing:

1. On each node, set the number of alignment threads with ``-t/--threads``,
   filtering or matrix loader workers with ``--nproc``, and ICE threads with
   ``--ice-threads``. Coarsening uses one process.
2. On an HPC cluster, you can submit the same ``runHiC pileup`` command on
   multiple nodes so separate jobs can work on different FASTQ chunks.

Example
-------
To run the :doc:`quickstart` example across nodes, create a shared working
directory and copy the metadata file into it::

    $ mkdir run-on-cluster
    $ cp datasets.tsv ./run-on-cluster
    $ cd run-on-cluster

Prepare the chromosome sizes and aligner index first. Submit the following
command from the shared directory on each node (for example, in a Slurm
job script). Use the same inputs and parameters, and give each job its own
log and temporary directory::

    $ runHiC pileup -p ../../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 20 --include-readid --drop-seq --chunkSize 1500000 --nproc 20 --ice-threads 20 -C ../../data/hg38/hg38.chrom.sizes -i ../../data/hg38/hg38.fa --tmpdir ".runHiC-${SLURM_JOB_ID}" --logFile "runHiC-${SLURM_JOB_ID}.log"

After the jobs finish, check that the expected matrices were produced. If
you want to run the stages separately, submit multiple ``mapping`` jobs in
the shared directory, then run ``filtering`` and ``binning`` once all mapping
jobs have finished.
