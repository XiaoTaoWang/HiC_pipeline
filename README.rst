Introduction
============
runHiC is a easy-to-use Hi-C processing software based on hiclib (https://bitbucket.org/mirnylab/hiclib) .
Different from hiclib, which was born for flexibility, runHiC is a customized pipeline, and can be
run from command line directly.

Links
=====
- `Repository <https://github.com/XiaoTaoWang/HiC_pipeline>`_
- `PyPI <https://pypi.python.org/pypi/runHiC>`_

Installation
============
Python Packages
---------------
We recommend using `conda <http://conda.pydata.org/miniconda.html>`_, an excellent Python package and
environment manager.

Open a terminal and type::

    $ conda install numpy numexpr scipy matplotlib cython biopython h5py pysam pip

Install bx-python and joblib using pip::

    $ pip install joblib bx-python

Install mirnylib and hiclib from source code:

Download `mirnylib <https://bitbucket.org/mirnylab/mirnylib>`_ and `hiclib <https://bitbucket.org/mirnylab/hiclib>`_,
and run install_linux.py contained in the unpacked folder, respectively.

Other dependencies
------------------
Install samtools:

Download `samtools <http://sourceforge.net/projects/samtools/files/>`_, unpack it, change to the extracted
directory::

    $ make

Make *samtools* accessible to the system. (Via environment variable *PATH*)

Install Bowtie2:

Download the `source code <http://sourceforge.net/projects/bowtie-bio/files/bowtie2/>`_, unzip it and
add the path to the extracted directory to *PATH*.

runHiC
------

Use easy_install::

    $ easy_install runHiC

Design Concepts
===============
runHiC is able to perform the entire analysis from sequencing data to corrected HeatMaps.

runHiC separate the whole process into 4 stages and you can begin and end at any stage using certain
subcommands.

7 subcommands are available:

- *mapping*: Iteratively map pair-end sequencing reads to a supplied genome
- *filtering*: Remove noises at the level of aligned read pairs and restriction fragments
- *binning*: Bin filtered reads at certain resolution (original Heat Maps are generated)
- *correcting*: Perform iterative corrections on the original Heat Maps
- *pileup*: Streamline all 4 subcommands above from *mapping* to *correcting*.
- *tosparse*: Convert intra-chromosomal contact matrices to sparse ones.

Preparation
===========
Before running this program, you need to carry out several other things to improve performance:

**Re-organize your directory arrangements**

Although not required, we recommend creating a data root directory separate from the working
directory.

**Place genome and sequencing data under the data root directory**

Genome sequences should be stored chromosome by chromosome in FASTA format under a subfolder named
after corresponding genome name. Sequencing read-pairs should be stored in SRA or FASTQ format under
another subfolder (any valid name).

**Construct a metadata file describing your sequencing data under the working directory**

Four columns are required: prefix of SRA file name, cell line name, biological replicate label, and
restriction enzyme name. An example file is distributed along with this software, please check it.

Usage
=====
Open a terminal, type ``runHiC -h`` and ``runHiC <subcommand> -h`` for help information.
