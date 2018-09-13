Quick Start
***********
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
current working directory to the sub-folder *demo*::

    $ cd demo
    $ mkdir data
    $ ls -lh

    total 16K
    drwxr-x--- 2 xtwang CPeng 4.0K Oct 29 10:03 data
    drwxr-x--- 2 xtwang CPeng 4.0K Oct 29 10:01 workspace

During this tutorial, all input data including Hi-C raw sequencing data and
the reference genome data will be placed under the *data* sub-folder, and
*runHiC* will be run under the *workspace* sub-folder.

Download the example Hi-C data set from a human cell line GM06990::

    $ cd data
    $ mkdir HiC-SRA
    $ cd HiC-SRA
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX011/SRX011608/SRR027956/SRR027956.sra -O SRR027956.sra
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX011/SRX011610/SRR027958/SRR027958.sra -O SRR027958.sra
    $ ls -lh

    total 1.4G
    -rw-r----- 1 xtwang CPeng 623M Oct 29 10:07 SRR027956.sra
    -rw-r----- 1 xtwang CPeng 783M Oct 29 10:10 SRR027958.sra

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
    $ cd hg19
    $ wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*
    $ wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz

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

Finally, uncompress the .gz files to finish this section::

    $ gunzip *.gz
    $ cd ../..
	
Mapping
=======
The first processing stage of *runHiC* is conducted by the *mapping* subcommand,
which maps raw sequencing reads to the reference genome and assigns aligned
reads to the restriction fragments.

*runHiC* also records the read-level statistics at this stage for quality
assessment of your Hi-C data. (See `quality <http://xiaotaowang.github.io/HiC_pipeline/quality.html>`_)

Usage
-----
``runHiC mapping [options]``

Bowtie2 Path
------------
*runHiC* calls the *bowtie2* read alignment software during the *mapping* stage,
so you need to tell *runHiC* where *bowtie2* is installed in your system::

    $ which bowtie2

    ~/Tools/anaconda2/bin/bowtie2

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

    $ runHiC mapping -m datasets.tsv -p ../data -g hg19 -G gap.txt -f HiC-SRA -F SRA -b ~/Tools/anaconda2/bin/bowtie2 -t 20 --removeInters --logFile runHiC-mapping.log

For FASTQ and the compressed FASTQ format, just replace "HiC-SRA" with "HiC-FASTQ" or "HiC-gzip",
and reset "-F" argument correspondingly::

    $ runHiC mapping -m datasets.tsv -p ../data -g hg19 -G gap.txt -f HiC-gzip -F FASTQ -b ~/Tools/anaconda2/bin/bowtie2 -t 20 --removeInters --logFile runHiC-mapping.log

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

  Path to the bowtie2 executable, see `Bowtie2 Path <http://xiaotaowang.github.io/HiC_pipeline/quickstart.html#bowtie2-path>`_.

- ``-t/--threads``

  Number of the bowtie2 threads.

- ``--removeInters``

  Remove intermediate results if specified.

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
>>> lib = h5dict.h5dict('./hdf5-hg19/SRR027956/SRR027956.hdf5','r')
>>> lib.keys()
[u'chrms1', u'chrms2', u'cuts1', u'cuts2', u'downrsites1', u'downrsites2', u'misc', u'rfragIdxs1', u'rfragIdxs2', u'rsites1', u'rsites2', u'strands1', u'strands2', u'uprsites1', u'uprsites2']

Filtering
=========
The *filtering* subcommand of *runHiC* is designed to perform some basic filtering on
the aligned read pairs: [1]_

1. Remove the read pair that maps to the same restriction fragment.
2. Remove redundant PCR artifacts.

During the filtering process, *runHiC* also records the fragment-level and the
contact-level statistics for quality assessment of your Hi-C data. (See `quality <http://xiaotaowang.github.io/HiC_pipeline/quality.html>`_)

Data from the same biological replicate (or optionally all replicates of the same cell
line) are merged in this processing stage.

Here's the command you should type in the terminal::

    $ runHiC filtering -m datasets.tsv --HDF5 hdf5-hg19 --libSize 500 --duplicates -l 2

- ``-m/--metadata``

  The metadata data file name.

- ``--HDF5``

  Path to the root folder of HDF5 files generated in the *mapping* stage.

- ``--libSize``

  Maximum length of molecules in your Hi-C library.
  
- ``--duplicates``

  Remove redundant PCR artifacts if specified.

- ``-l/--level``

  Merging level. 1: Merge data from the same biological replicate; 2: Merge data from
  all replicates of the same cell line.

That will create a new sub-folder named *filtered-hg19* in which the filtered data
(stored in the HDF5 format, suffixed with .hdf5) reside. (See `Access Data from HDF5 <http://xiaotaowang.github.io/HiC_pipeline/quickstart.html#access-data-from-hdf5>`_)

Binning
=======
This processing stage is conducted by the *binning* subcommand.

By partitioning the linear genome into fixed-size bins (intervals), the number of contacts
(the filtered read pairs) can be counted between any two bins (reads are assigned to bins
by the centers of their fragments), which results in a "contact" matrix.

*runHiC* supports two modes of binning: **wholeGenome** builds a single genome-versus-genome
contact matrix, and **byChromosome** builds chromosome-versus-chromosome contact matrices.

Type in the command below to carry on our analyzing on the example data::

    $ runHiC binning -f filtered-hg19 -M byChromosome -R 2000000

- ``-f/--filteredDir``

  Path to the hdf5 files generated during the *filtering* stage. Wild cards are allowed.
  If a path points to a folder, the *binning* procedure will be performed on each hdf5
  file under that folder.

- ``-M/--mode``

  Mode label for building contact matrices. *wholeGenome* or *byChromosome*.

- ``-R/--resolution``

  Resolution of the generated contact matrices. Unit: bp

After this command, a new sub-folder named *Raw-hg19* will be created under current
working directory. The original contact matrices are also stored in hdf5 files suffixed
with ".hm". (See `Access Data from HDF5 <http://xiaotaowang.github.io/HiC_pipeline/quickstart.html#access-data-from-hdf5>`_)

Correcting
==========
Hi-C data can contain many different biases, some from known origin and others from an
unknown origin. There are two general approaches to Hi-C bias correction: explicit-factor
methods [2]_ and matrix balancing methods [1]_. *runHiC* (hiclib) uses the matrix balancing
algorithm called Sinkhorn–Knopp for bias correction.

To correct our sample Hi-C data, type in the command below::

    $ runHiC correcting -H Raw-hg19

- ``-H/--Heatmap``

  Path to the original contact matrices generated during the last stage. Wild cards are
  allowed. If a path points to a folder, bias correction will be performed on each ".hm"
  file under that folder.

After that, a new sub-folder named *Corrected-hg19* with corrected contact matrices will
be created under current working directory.

tosparse
========
Contact matrices generated by *binning* and *correcting* subcommand are represented as
dense matrices, which are memory-consuming for large genomes (such as human and mouse)
and high resolutions (40 Kb) when used in other calculations (such as TAD and loop identifications).
This subcommand converts intra-chromosomal contact matrices into sparse ones.

To convert the dense matrices generated just now into well-known CSR (Compressed Sparse
Row) matrices, type in the command below::

    $ runHiC tosparse -H Corrected-hg19 --csr
	
-- ``-H/--cHeatMap``

  Path to the dense matrix files from *binning* or *correcting*. Wild cards are allowed.
  If a folder name is provided, conversion will be performed on each ".hm" file under
  that folder.

-- ``--csr``

  If specified, dense matrices are converted into `CSR (Compressed Sparse Row) <https://docs.scipy.org/doc/scipy/reference/sparse.html>`_ matrices,
  a `customized numpy structured array <http://xiaotaowang.github.io/TADLib/hitad.html#transform-txt-into-npz>`_ is applied otherwise.

Sparse matrices are organized in `npz <https://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html>`_ format
using chromosome labels as the keys. Here, you can find them under the "Corrected-hg19" folder.

Pileup
======
*runHiC* also provides a handy subcommand called "pileup" by which you can perform all
processing steps above with single-line command::

    $ runHiC pileup -p ../data -g hg19 -f HiC-SRA -F SRA -b ~/Tools/anaconda2/bin/bowtie2 -t 20 --removeInters -M byChromosome -R 2000000



Reference
=========
.. [1] Imakaev M, Fudenberg G, McCord RP et al. Iterative correction of Hi-C data
       reveals hallmarks ofchromosome organization. Nat Methods, 2012, 9(10): 999-1003.

.. [2] Yaffe E, Tanay A. Probabilistic modeling of Hi-C contact maps eliminates
       systematic biases to characterize global chromosomal architecture. Nat Genet,
       2011, 43(11): 1059-65.

