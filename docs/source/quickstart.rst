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
Suppose you are still in the *runHiC* distribution root folder::

    $ ls -lh

?

Create a new sub-folder named *demo* and change you current working directory
to it::

    $ mkdir demo
	$ cd demo
	$ pwd

?

Create a sub-folder named *data*. Both Hi-C raw sequencing read data and the
reference genome data will be placed under it::

    $ mkdir data
	$ cd data
	$ pwd

?

Download the example Hi-C data set from a human cell line GM06990::

    $ mkdir HiC-SRA
	$ cd HiC-SRA
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX011/SRX011608/SRR027956/SRR027956.sra -O SRR027956.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX011/SRX011610/SRR027958/SRR027958.sra -O SRR027958.sra
	$ ls -lh

?

Reads in SRR027956.sra and SRR027958.sra are sequenced from two biological replicates,
respectively.

To demonstrate all read formats supported by runHiC, let's first dump reads
by *fastq-dump* and then compress the FASTQ files using gzip::

    $ cd ..
	$ mkdir HiC-FASTQ
	$ cd HiC-FASTQ
	$ for i in ../HiC-SRA/*.sra; do fastq-dump --split-3 $i; done
	$ for i in ./*.fastq; do gzip -c $i >


	