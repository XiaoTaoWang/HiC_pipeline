runHiC
******
runHiC is an easy-to-use command-line tool for Hi-C data processing.

Since versio 0.8.5, runHiC has changed the default aligner to `chromap <https://github.com/haowenz/chromap>`_,
which is comparable to `bwa-mem <https://github.com/lh3/bwa>`_ in alignment accuracy, but runs over 10 times faster.

Since version 0.8.1, runHiC can be used directly on `Arima HiC <https://arimagenomics.com>`_ data
by setting the enzyme name to *Arima*.

Since version 0.8.0, runHiC has changed its default data container/format from HDF5 to
`Pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_ and
`Cooler <https://github.com/mirnylab/cooler>`_. 

Design Concepts
===============
runHiC is designed to process Hi-C data from raw sequencing reads(.sra, .fastq, .fastq.gz) to the ICE-corrected
contact matrices. It currently contains 5 subcommand:

+------------+-------------------------------------------------------------------------------------------------------------------+
| mapping    | Map raw pair-end sequencing data to a supplied genome. Support three read aligners: chromap, bwa and minimap2.    |
+------------+-------------------------------------------------------------------------------------------------------------------+
| filtering  | Perform read-level and fragment-level noise removing                                                              |
+------------+-------------------------------------------------------------------------------------------------------------------+
| binning    | 1.Generate contact matirx; 2. Perform ICE/matrix-balancing normalization                                          |
+------------+-------------------------------------------------------------------------------------------------------------------+
| pileup     | Perform entire processing from *mapping* to *binning*                                                             |
+------------+-------------------------------------------------------------------------------------------------------------------+
| quality    | Assess the quality of your Hi-C data                                                                              |
+------------+-------------------------------------------------------------------------------------------------------------------+


User Guide
==========

.. toctree::
   :maxdepth: 3

   install
   quickstart
   parallel
   quality
   changelog

Citation
========
Xiaotao Wang. (2016). runHiC: A user-friendly Hi-C data processing software based on hiclib. Zenodo.
`10.5281/zenodo.55324 <http://dx.doi.org/10.5281/zenodo.55324>`_

