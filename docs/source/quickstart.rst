Quick Start
***********
This tutorial covers the basic usage of 4 subcommands (*mapping*, *filtering*,
*binning*, and *pileup*) provided by runHiC. We will first download an example
Hi-C dataset and corresponding reference genome data. Then we will process the
Hi-C data step by step from raw sequencing reads (.sra, .fastq, and .fastq.gz)
to the ICE-corrected contact matrices. Lastly, we will demonstrate how to streamline
the processing pipeline by using the *pileup* subcommand.

Data Preparation
================
First of all, let's make a temporary blank folder somewhere in your system, change your
working directory to it, and make two sub-folders named *data* and *workspace* within it::

    $ mkdir data
    $ mkdir workspace
    $ ls -lh

    total 0
    drwxr-xr-x  2 xtwang  staff    64B Sep 16 11:24 data
    drwxr-xr-x  2 xtwang  staff    64B Sep 16 11:25 workspace

During this tutorial, all input data including the raw sequencing data and
the reference genome data will be placed under the *data* sub-folder, and
*runHiC* will be run under the *workspace* sub-folder.

Then download an example Hi-C dataset using the *prefetch* command of the SRA toolkit::

    $ cd data
    $ mkdir HiC-SRA
    $ cd HiC-SRA
    $ prefetch -o SRR027956.sra SRR027956
    $ prefetch -o SRR027958.sra SRR027958
    $ ls -lh

    total 1.4G
    -rw-r--r-- 1 xtwang staff 623M Sep 16 11:49 SRR027956.sra
    -rw-r--r-- 1 xtwang staff 783M Sep 16 11:53 SRR027958.sra

*runHiC* currently supports three read formats: SRA (Sequence Read Archive), FASTQ,
and compressed FASTQ (.fastq.gz). To demonstrate this, let's first dump reads
using the *fastq-dump* command, and then compress the FASTQ files using *gzip*::

    $ for i in ./*.sra; do fastq-dump --split-files $i; done
    $ for i in ./*.fastq; do gzip -c $i > `basename $i`.gz; done
    $ cd ..
    $ mkdir HiC-FASTQ
    $ mkdir HiC-gzip
    $ mv ./HiC-SRA/*.fastq ./HiC-FASTQ
    $ mv ./HiC-SRA/*.gz ./HiC-gzip

.. note::
   If your FASTQ files are suffixed with "R1.fastq" and "R2.fastq", please make
   sure rename them as "_1.fastq" and "_2.fastq" before you run runHiC, as if
   they are dumped from an sra file.
	
Then download the reference genome (hg38) data from UCSC::

    $ mkdir hg38
    $ cd hg38
    $ wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/*

Note that above commands can be modified to download any other genomes available
in UCSC, by replacing "hg38" with the desired reference genome release name.

Let's include chromosomes that are completely assembled only. To do so, open a python
interpreter and follow the commands below:

>>> import os, glob
>>> labels = list(map(str,range(1,23))) + ['X','Y','M']
>>> pool = ['chr{0}.fa.gz'.format(i) for i in labels]
>>> for c in glob.glob('*.fa.gz'):
...     if not c in pool:
...         os.remove(c)
>>> exit()

Finally, uncompress the .gz files and merge all chromosomes into hg38.fa::

    $ gunzip *.gz
    $ cat *.fa > hg38.fa
    $ cd ../..
	
Mapping
=======
The first step of *runHiC* is conducted by the *mapping* subcommand,
which maps raw sequencing reads to the reference genome and parses the resulted
SAM/BAM alignments into `.pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_.

Usage
-----
``runHiC mapping [options]``

Create the Meta Data File
-------------------------
Before running *runHiC*, another thing you need to do is to create a TXT file
named "datasets.tsv" under the *workspace* sub-folder::

    $ cd workspace
    $ vim datasets.tsv

And fill in the following content::
    
    SRR027956 GM06990 R1 HindIII
    SRR027958 GM06990 R2 HindIII
	
Here, "datasets.tsv" is a meta data file describing your Hi-C data, which should
contain 4 columns. In order, they are: the prefix of the SRA file name (in the
case of the FASTQ read format, it should be the leading part of the file names
apart from the "_1.fastq" or "_2.fastq" substring), cell line name, biological
replicate label, and the restriction enzyme name. Note that for Arima Hi-C, you
can set the enzyme name to *Arima*; for experiments that do not use restriction
enzymes for DNA fragmentation, you can set the enzyme name arbitrarily for your
record. For example, for Micro-C, you can set it to *MNase*; for ChIA-PET, you
can set it to *sonication*.

runHiC Command
---------------
Now type and execute the command below::

    $ runHiC mapping -p ../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 10 --include-readid --drop-seq --chunkSize 1500000 --logFile runHiC-mapping.log

For FASTQ and the compressed FASTQ format, replace "HiC-SRA" with "HiC-FASTQ"
or "HiC-gzip", and reset the "-F" argument accordingly::

    $ runHiC mapping -p ../data/ -g hg38 -f HiC-gzip -F FASTQ -A bwa-mem -t 10 --include-readid --drop-seq --chunkSize 1500000 --logFile runHiC-mapping.log

Two sub-folders named *alignments-hg38* and *pairs-hg38* will be created under current
working directory (*workspace*). During this process:

1. Read pairs will be mapped to the *hg38* reference genome with ``bwa mem``, and the
   alignment results will be reported in BAM format and placed under *alignments-hg38*.
2. BAM files will be parsed into `.pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_
   using `pairtools <https://github.com/mirnylab/pairtools>`_ under *pairs-hg38*.

*runHiC* currently supports three read aligners, *bwa-mem*, *chromap*, and *minimap2*.
You can switch it by the ``-A/--aligner`` argument.

During the alignment parsing, *runHiC* detects ligation junctions, marks various situations
(Unmapped, Multimapped, Multiple ligations-Walks, and valid Single ligations), and sorts
pairs for further analysis. In this example, .pairsam.gz files under *pairs-hg38* are
valid .pairs files defined by the `4DN <https://www.4dnucleome.org>`_ consortium. By default,
it will only contain 7 columns: chr1, pos1, chr2, pos2, strand1, strand2, and pair_type;
if you add ``--include-readid`` on the command, you will get an additional "readID" column;
if you specify ``--include-sam``, two extra columns "sam1" and "sam2" will be added to store
the original alignments; if you add ``--drop-seq``, SEQ and QUAL will be removed from the sam
fields to save the disk space.

Filtering
=========
The *filtering* subcommand of *runHiC* is designed to perform basic filtering procedures on
the aligned read pairs. These filtering procedures include:

1. Remove redundant PCR artifacts.
2. Remove the read pair that maps to the same restriction fragment (since version 0.8.5, runHiC
   only performs this filtering if you specify ``--add-frag`` when you run ``runHiC mapping``).

During the filtering process, *runHiC* also records read-level, fragment-level and the
contact-level statistics for quality assessment of your Hi-C data.
(See `quality <http://xiaotaowang.github.io/HiC_pipeline/quality.html>`_)

Here's the command you should type in the terminal::

    $ runHiC filtering --pairFolder pairs-hg38/ --logFile runHiC-filtering.log --nproc 10

That will create a new sub-folder named *filtered-hg38*. Please find the final valid
contact pairs in .pairs.gz files. If you specify ``--include-sam`` when you run
``runHiC mapping``, it will also output a .bam file accompanying each .pairs.gz file
to store alignments that passed all filtering criteria.

Binning
=======
In this step, an .mcool file will be produced under the *coolers-hg38* sub-folder for each
.pairs.gz file using `cooler <https://cooler.readthedocs.io/en/latest/>`_. The mcool format
is the official Hi-C data format for the `4DN consortium <https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline>`_
and can be visualized using `HiGlass <https://docs.higlass.io/>`_::

    $ runHiC binning -f filtered-hg38/ --logFile runHiC-binning.log --nproc 10

Pileup
======
*runHiC* also provides a handy subcommand called "pileup" by which you can perform all
processing steps above using the single-line command below::

    $ runHiC pileup -p ../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 10 --include-readid --drop-seq --chunkSize 1500000 --logFile runHiC.log