Introduction
============
runHiC supports parallel processing of large Hi-C datasets on both SMP (systems
with multiple processors or cores) and PBS-based clusters (computers connected
via network).

All you need to know is another 3 command-line arguments (refer to `quickstart?`
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

???

Now submit your PBS script by using the ``qsub`` command::

???

That's it! Just a single-line command and no tedious setup required.

Performance
===========
The table below shows the performance test of runHiC with low-depth and high-depth sequencing
data. (Running time means the wall time runHiC took from raw read files to corrected sparse
contact matrices)

+----------------+-------------------+------------------+-----------------------------+
| Cell Line      | Number of reads   | Number of nodes  | Running time (hr: min: sec) |
+================+===================+==================+=============================+



