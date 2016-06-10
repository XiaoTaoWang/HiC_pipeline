runHiC
******
.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.55324.svg
   :target: http://dx.doi.org/10.5281/zenodo.55324

Introduction
============
runHiC is a easy-to-use Hi-C processing software based on hiclib (https://bitbucket.org/mirnylab/hiclib)
Different from hiclib, which was born for flexibility, runHiC is a customized pipeline, and can be
run from command line directly.

Links
=====
- `Repository <https://github.com/XiaoTaoWang/HiC_pipeline>`_
- `PyPI <https://pypi.python.org/pypi/runHiC>`_

Installation
============
Please check the file "INSTALL.rst" in the distribution.

Design Concepts
===============
runHiC is able to perform the entire analysis from sequencing data to corrected contact matrices. It
separates the whole process into 4 stages(*mapping*, *filtering*, *binning*, *correcting*). You can
begin and end at any stage using certain subcommands.

7 subcommands are available:

+------------+------------------------------------------------------------------------------+
| mapping    | Iteratively map pair-end sequencing reads to a supplied genome               |
+------------+------------------------------------------------------------------------------+
| filtering  | Remove noises at the level of aligned read pairs and restriction fragments   |
+------------+------------------------------------------------------------------------------+
| binning    | Generate original contact matrices                                           |
+------------+------------------------------------------------------------------------------+
| correcting | Perform iterative corrections on original contact matrices                   |
+------------+------------------------------------------------------------------------------+
| tosparse   | Convert intra-chromosomal contact matrices to sparse ones                    |
+------------+------------------------------------------------------------------------------+
| pileup     | Streamline all stages from *mapping* to *correcting*                         |
+------------+------------------------------------------------------------------------------+
| quality    | Assess the quality of your experiments                                       |
+------------+------------------------------------------------------------------------------+
| visualize  | Plot the heatmap for given interval                                          |
+------------+------------------------------------------------------------------------------+

Preparation
===========
Please refer to the **Sample** folder distributed with our source code.

Directory Rearrangements
````````````````````````
Although not required, I recommend creating a data root directory separate from the working
directory.

Data Placement
``````````````
Both genome and sequencing data should be placed under the data root directory.

Genome sequences should be stored chromosome by chromosome in FASTA format under a subfolder(named
after corresponding genome name).

Sequencing read-pairs should be stored in SRA or FASTQ format under another subfolder(any valid name).

Meta Data
`````````
Construct a meta data file describing your sequencing data under the working directory

Four columns are required: prefix of SRA file name, cell line name, biological replicate label, and
restriction enzyme name. An example file(Sample/working/datasets.tsv) is distributed along with this
software, please check it.

Usage
=====
Open a terminal, type ``runHiC -h`` and ``runHiC <subcommand> -h`` for help information.
