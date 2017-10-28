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

  The maximum number of task processes to launch on a single machine.

- ``--chunkSize``

  

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
The table below shows the performance test of runHiC with low-thoughput and high-thoughput
