Introduction
============
runHiC supports parallel processing of large Hi-C datasets on both SMP (systems
with multiple processors or cores) and PBS-based clusters (computers connected
via network).

All you need to know is another 3 command-line arguments (refer to `quickstart <http://xiaotaowang.github.io/HiC_pipeline/quickstart.html>`_
for basic usage):

- ``-r/--running-mode``

  Running mode of the program. "local" (Run on a local computer) or "pbs" (Run on
  a PBS-based cluster).

- ``--nworker``

  The maximum number of task processes to launch on a single machine (which means
  the local computer for "local" mode and a compute node for "pbs" mode).

- ``--chunkSize``

  To make parallelized computation possible, our pipeline splits the original read file into
  chunks and maps them to the reference genome separately. This parameter specifies the size
  of each chunk. By default, no split is performed.

  

Example
=======
To show how to launch runHiC on our example data (see `Data Preparation`) across
multiple nodes of a cluster, let's create a sub-folder and copy the meta data file
into it::

    $ mkdir run-on-cluster
    $ cp datasets.tsv ./run-on-cluster

Then change to the sub-folder and write (vi) a simple PBS script. My job script is
shown below::

    #PBS -N runHiC
    #PBS -l nodes=4:ppn=20
    #PBS -l walltime=50:00:00
    #PBS -q batch
    #PBS -V
    #PBS -S /bin/bash

    cd $PBS_O_WORKDIR

    runHiC pileup -r pbs --nworker 2 -p ../../data -g hg19 -f HiC-gzip -F FASTQ -b ~/Tools/anaconda2/bin/bowtie2 -t 20 --chunkSize 2000000 -M byChromosome -R 2000000 --logFile runHiC.log

Now submit your PBS script by using the ``qsub`` command::

    $ qsub pileup.pbs

That's it! Just a single-line command and no tedious setup required.

Performance
===========
The table below shows the performance test of runHiC with low-depth and high-depth sequencing
data. (Running time means the wall time runHiC took from raw read files to corrected sparse
contact matrices)

+----------------+-------------------+------------------+-----------------------------+
| Cell Line      | Number of reads   | Number of nodes  | Running time (hr: min: sec) |
+================+===================+==================+=============================+
| GM06990        | 14,332,993        |        1         |      0:56:01                |
+----------------+-------------------+------------------+-----------------------------+
| GM06990        | 14,332,993        |        4         |      0:16:58                |
+----------------+-------------------+------------------+-----------------------------+
| MCF7           | 944,178,846       |        10        |     11:14:24                |
+----------------+-------------------+------------------+-----------------------------+
