Introduction
============
This tutorial covers the basic use of 6 subcommands (*mapping*, *filtering*,
*binning*, *correcting*, *tosparse*, and *pileup*) provided by runHiC. We
will first download an example Hi-C data set and corresponding reference
genome data. Then we will use runHiC to process the Hi-C data step by step
from raw sequencing reads (.sra, .fastq, .fastq.gz) to the corrected
contact matrices. Lastly, we will demonstrate how to streamline the processing
pipeline by *pileup* subcommand.

Data Preparation
================
Suppose you are still in the *runHiC* distribution root folder, change your
current working directory to the sub-folder **demo**::

    $ cd demo
	$ ls -lh

?

During this tutorial, all input data including Hi-C raw sequencing data and
the reference genome data will be placed under the **data** sub-folder, and
*runHiC* will be run under the *workspace* sub-folder.

Download the example Hi-C data set from a human cell line GM06990::

    $ cd data
    $ mkdir HiC-SRA
	$ cd HiC-SRA
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX011/SRX011608/SRR027956/SRR027956.sra -O SRR027956.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX011/SRX011610/SRR027958/SRR027958.sra -O SRR027958.sra
	$ ls -lh

?

Reads in SRR027956.sra and SRR027958.sra are sequenced from two biological replicates,
respectively.

To demonstrate all read formats supported by *runHiC*, let's first dump reads
with *fastq-dump* and then compress the FASTQ files with gzip::

	$ for i in ./*.sra; do fastq-dump --split-3 $i; done
	$ for i in ./*.fastq; do gzip -c $i > `basename $i`.gz; done

*runHiC* currently supports 3 sequencing read format: FASTQ, compressed FASTQ with gzip,
and SRA(Sequence Read Archive)::

    $ cd ..
	$ mkdir HiC-FASTQ
	$ mkdir HiC-gzip
	$ mv ./HiC-SRA/*.fastq ./HiC-FASTQ
	$ mv ./HiC-SRA/*.gz ./HiC-gzip
	
Download the reference genome (hg19) data from UCSC::

    $ mkdir hg19
	$ cd  hg19
    $ rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ .
	$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz

Above commands can be modified to download the data for any other genomes
available in UCSC by replacing "hg19" with the desired reference genome
release name.

We will only include sequences that have been completely assembled. Open
a python interpreter and follow the commands below:

>>> import os, glob
>>> labels = map(str,range(1,23)) + ['X','Y','M']
>>> pool = ['chr{0}.fa.gz'.format(i) for i in labels]
>>> for c in glob.glob('*.fa.gz'):
...     if not c in pool:
...         os.remove(c)
>>> exit()

Finally, uncompress the gz files to finish this section::

	$ gunzip *.gz
	$ cd ../..
	
Mapping
=======
The first processing stage of *runHiC* is conducted by the *mapping* subcommand,
which maps raw sequencing reads to the reference genome and assigns aligned
reads to the restriction fragments.

*runHiC* also records the read-level statistics at this stage for quality
assessment of your Hi-C data. (See ?)

Usage
-----
``runHiC mapping [options]``

.. _locate-the-bowtie2

Bowtie2 Path
------------
*runHiC* calls the *bowtie2* read alignment software during the *mapping* stage,
so you need to tell *runHiC* where *bowtie2* is installed in your system::

    $ which bowtie2

?

Write the Meta Data File
------------------------
Another thing you need to do is to prepare a meta data file describing your Hi-C
data under the *workspace* sub-folder::

    $ cd workspace
	$ cat datasets.tsv
	
The meta data file should contain 4 columns: prefix of the SRA file name (in the
case of the FASTQ read format, it should be the leading part of the file name
apart from the "_1.fastq" or "_2.fastq" substring), cell line name, biological
replicate label, and the restriction enzyme name::

    SRR027956 GM06990 R1 HindIII
    SRR027958 GM06990 R2 HindIII

runHiC Command
---------------
Now type in the command below::

    $ ?

- ``-m/--metadata``

  The metadata data file name.

- ``-p/--dataFolder``

  Path to the root folder containing both Hi-C sequencing data and the reference
  genome data.

- ``-g/--genomeName``

  Name of the reference genome. (Or name of the folder containing the reference
  genome data)

- ``-G/--gapFile``

  Name of the decompressed gap file downloaded from UCSC. If runHiC fails to find
  it, a dummy one will be generated in the specified genome folder (see ``-g/--genomeName``).

- ``-f/--fastqDir``

  Name of the folder containing the Hi-C raw sequencing data.

- ``-F/Format``

  Format of the sequencing data. SRA or FASTQ.

- ``-b/--bowtiePath``

  Path to the bowtie2 executable, see :ref:`locate-the-bowtie2`.

- ``-t/--threads``

  Number of the bowtie2 threads.

- ``--removeInters``

  Whether to remove intermediate results.

- ``--logFile``

  Log file name.

During the execution of ``runHiC mapping``, two new sub-folders named *bams-hg19* and
*hdf5-hg19* are created under current working directory (*workspace*). The read pairs
are mapped to the *hg19* reference genome in an iterative way with *bowtie2*. [1]_
The alignment results are stored in the BAM format and placed under *bams-hg19*. Then
BAM files of corresponding read pairs are parsed together and outputed into HDF5 files
(suffixed with .hdf5) under *hdf5-hg19*.

runHiC uses a rotating file for logging. According to our settings, when the size of
"runHiC.log" gets about 100K, it is closed and renamed to "runHiC.log.1". At the same
time, a new file "runHiC.log" is silently opened for output. In a word, the system saves
old log files by appending the extensions ".1", ".2" etc., and the current log is always
written to "runHiC.log".

Access Data from HDF5
---------------------
You can extract data from HDF5 files via *mirnylib*:

>>> from mirnylib import h5dict
>>> ???

Filtering
=========
The *filtering* subcommand of *runHiC* is designed to perform some basic filtering on
the aligned read pairs: [1]_

1. Remove the read pair that maps to the same restriction fragment.
2. Remove redundant PCR artifacts.

During the filtering process, *runHiC* also records the fragment-level and the
contact-level statistics for quality assessment of your Hi-C data. (See ?)

Data from the same biological replicate (or optionally all replicates of the same cell
line) will be merged in this processing stage.




Reference
=========
.. [1] Imakaev M, Fudenberg G, McCord RP et al. Iterative correction of Hi-C data
       reveals hallmarks ofchromosome organization. Nat Methods, 2012, 9: 999-1003.

