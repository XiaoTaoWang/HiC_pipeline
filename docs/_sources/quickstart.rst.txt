Quick Start
***********
This tutorial covers the basic use of 4 subcommands (*mapping*, *filtering*,
*binning*, and *pileup*) provided by runHiC. We will first download an example
Hi-C data set and corresponding reference genome data. Then we will use runHiC
to process the Hi-C data step by step from raw sequencing reads (.sra, .fastq, .fastq.gz)
to the corrected contact matrices. Lastly, we will demonstrate how to streamline
the processing pipeline by *pileup* subcommand.

Data Preparation
================
First of all make a temporary blank folder somewhere on your system, and change your
working directory to it::

    $ mkdir data
    $ mkdir workspace
    $ ls -lh

    total 0
    drwxr-xr-x  2 xtwang  staff    64B Sep 16 11:24 data
    drwxr-xr-x  2 xtwang  staff    64B Sep 16 11:25 workspace

During this tutorial, all input data including Hi-C raw sequencing data and
the reference genome data will be placed under the *data* sub-folder, and
*runHiC* will be run under the *workspace* sub-folder.

Download the example Hi-C data set using *prefetch* of SRA toolkit::

    $ cd data
    $ mkdir HiC-SRA
    $ cd HiC-SRA
    $ prefetch -O . SRR027956 
    $ prefetch -O . SRR027958
    $ ls -lh

    total 1.4G
    -rw-r--r-- 1 xtwang staff 623M Sep 16 11:49 SRR027956.sra
    -rw-r--r-- 1 xtwang staff 783M Sep 16 11:53 SRR027958.sra

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
	
Download the reference genome (hg38) data from UCSC::

    $ mkdir hg38
    $ cd hg38
    $ wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/*

Above commands can be modified to download any other genomes available in UCSC
by replacing "hg38" with the desired reference genome release name.

We will only include sequences that have been completely assembled. Open
a python interpreter and follow the commands below:

>>> import os, glob
>>> labels = list(map(str,range(1,23))) + ['X','Y']
>>> pool = ['chr{0}.fa.gz'.format(i) for i in labels]
>>> for c in glob.glob('*.fa.gz'):
...     if not c in pool:
...         os.remove(c)
>>> exit()

Finally uncompress the .gz files and merge all chromosomes into hg38.fa::

    $ gunzip *.gz
    $ cat *.fa > hg38.fa
    $ cd ../..
	
Mapping
=======
The first processing stage of *runHiC* is conducted by the *mapping* subcommand,
which just maps raw sequencing reads to the reference genome and parses the resulted
SAM/BAM alignments into `.pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_.

Usage
-----
``runHiC mapping [options]``

Write the Meta Data File
------------------------
Another thing you need to do is to prepare a meta data file describing your Hi-C
data under the *workspace* sub-folder::

    $ cd workspace
    $ vim datasets.tsv

Create a TXT file called "datasets.tsv" by ``vim`` and fill in the following content::
    
    SRR027956 GM06990 R1 HindIII
    SRR027958 GM06990 R2 HindIII
	
The meta data file should contain 4 columns: prefix of the SRA file name (in the
case of the FASTQ read format, it should be the leading part of the file name
apart from the "_1.fastq" or "_2.fastq" substring), cell line name, biological
replicate label, and the restriction enzyme name.

runHiC Command
---------------
Now type and execute the command below::

    $ runHiC mapping -p ../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 4 -O BAM --include-readid --include-sam --drop-seq --chunkSize 1500000 --logFile runHiC-mapping.log

For FASTQ and the compressed FASTQ format, just replace "HiC-SRA" with "HiC-FASTQ"
or "HiC-gzip", and reset "-F" argument correspondingly::

    $ runHiC mapping -p ../data/ -g hg38 -f HiC-gzip -F FASTQ -A bwa-mem -t 4 -O BAM --include-readid --include-sam --drop-seq --chunkSize 1500000 --logFile runHiC-mapping.log

Two sub-folders named *alignments-hg38* and *pairs-hg38* will be created under current
working directory (*workspace*):

1. Read pairs will be mapped to the *hg38* reference genome with ``bwa mem``, and the
   alignment results will be reported in BAM format and placed under *alignments-hg38*.
2. BAM files will be parsed into `.pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_
   using `pairtools <https://github.com/mirnylab/pairtools>`_ under *pairs-hg38*.

*runHiC* now supports two read aligner, *bwa* and *minimap2*. You can switch it by ``-A/--aligner``
argument. Generally, *bwa* performs better on short reads (<=100bp), and *minimap2* runs
faster and is similiarly accurate as *bwa* on longer reads (>100bp).

For alignment format, you can choose between SAM(TXT) and BAM(Binary) with ``-O/--outformat``
argument.

During alignment parsing, *runHiC* detects ligation junctions, marks various situations
(Unmapped, Multimapped, Multiple ligations-Walks, and valid Single ligations), and sort
pairs for further analysis. In this example, .pairsam.gz files under *pairs-hg38* are
valid .pairs files proposed by `4DN <https://www.4dnucleome.org>`_ group. By default,
it will only contain 7 columns: chr1, pos1, chr2, pos2, strand1, strand2 and pair_type;
if you add ``--include-readid`` on the command, you will get an additional "readID" column;
if you specify ``--include-sam``, two extra columns "sam1" and "sam2" will be added to store
the original alignments; if you use ``drop-seq``, SEQ and QUAL will be removed from the sam
fields.

*runHiC* uses a rotating file for logging. According to our settings, when the size of
"runHiC.log" gets about 100K, it is closed and renamed to "runHiC.log.1". At the same
time, a new file "runHiC.log" is silently opened for output. In a word, the system saves
old log files by appending the extensions ".1", ".2" etc., and the current log is always
written to "runHiC.log".

Filtering
=========
The *filtering* subcommand of *runHiC* is designed to perform some basic filtering on
the aligned read pairs:

1. Remove redundant PCR artifacts.
2. Remove the read pair that maps to the same restriction fragment.

During the filtering process, *runHiC* also records read-level, fragment-level and the
contact-level statistics for quality assessment of your Hi-C data.
(See `quality <http://xiaotaowang.github.io/HiC_pipeline/quality.html>`_)

Data from the same biological replicate and all replicates of the same cell line are merged
in this processing stage.

Here's the command you should type in the terminal::

    $ runHiC filtering --pairFolder pairs-hg38/ --genomepath ../data/hg38/hg38.fa

That will create a new sub-folder named *filtered-hg38*. Please find the final valid
contact pairs in *.pairs.gz files. If you specified ``--include-sam`` when you ran
``runHiC mapping``, it will also output a .bam file accompanying each .pairs.gz file
to store alignments that passed all filtering criteria.


Binning
=======
Since 0.8.0, *runHiC* has integrated *binning* with *correcting*, and you can perform
contact matrix building and ICE correcting [1]_ with a single subcommand *binning*::

    $ runHiC binning -f filtered-hg38/ -R 500000

After this command, a new sub-folder named *coolers-hg38* will be created under current
working directory. And a contact matrix at 500Kb resolution will be created and corrected
in `cooler <https://cooler.readthedocs.io/en/latest/>`_ format.

Pileup
======
*runHiC* also provides a handy subcommand called "pileup" by which you can perform all
processing steps above with single-line command::

    $ runHiC pileup -p ../data/ -g hg38 -f HiC-SRA -F SRA -A bwa-mem -t 4 -O BAM --include-readid --include-sam --drop-seq --chunkSize 1500000 -R 500000 --logFile runHiC.log



Reference
=========
.. [1] Imakaev M, Fudenberg G, McCord RP et al. Iterative correction of Hi-C data
       reveals hallmarks ofchromosome organization. Nat Methods, 2012, 9(10): 999-1003.

