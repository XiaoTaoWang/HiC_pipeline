runHiC
******
.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.55324.svg
   :target: http://dx.doi.org/10.5281/zenodo.55324
.. image:: https://static.pepy.tech/personalized-badge/runhic?period=total&units=international_system&left_color=black&right_color=orange&left_text=Downloads
   :target: https://pepy.tech/project/runhic

runHiC is an easy-to-use command-line tool for Hi-C data processing.

Since version 0.8.6, runHiC has supported data from all kinds of 3C-based experiments,
including Hi-C, Micro-C, HiChIP/PLAC-Seq, and ChIA-PET. For experiments that do not use
restriction enzymes for DNA fragmentation, you can set the enzyme name arbitrarily for your
record. For example, for Micro-C, you can set it to *MNase*; for ChIA-PET, you can set it to
*sonication*.

Since version 0.8.5, runHiC has changed the default aligner to `chromap <https://github.com/haowenz/chromap>`_,
which is comparable to `bwa-mem <https://github.com/lh3/bwa>`_ in alignment accuracy, but runs over 10 times faster.

Since version 0.8.1, runHiC can be used directly on `Arima HiC <https://arimagenomics.com>`_ data
by setting the enzyme name to *Arima*.

Since version 0.8.0, runHiC has changed its default data container/format from HDF5 to
`Pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_ and
`Cooler <https://github.com/mirnylab/cooler>`_. 

Design Concepts
===============
runHiC is designed to process Hi-C data from raw sequencing reads(.sra, .fastq, .fastq.gz) to the ICE-corrected
contact matrices. It currently contains 5 subcommands:

+------------+-------------------------------------------------------------------------------------------------------------------+
| mapping    | Map raw sequencing reads to a supplied genome. Support three read aligners: chromap, bwa and minimap2.            |
+------------+-------------------------------------------------------------------------------------------------------------------+
| filtering  | Perform read-level and fragment-level noise removing                                                              |
+------------+-------------------------------------------------------------------------------------------------------------------+
| binning    | 1.Generate contact matirx; 2. Perform ICE/matrix-balancing normalization                                          |
+------------+-------------------------------------------------------------------------------------------------------------------+
| pileup     | Perform the entire processing steps from *mapping* to *binning*                                                   |
+------------+-------------------------------------------------------------------------------------------------------------------+
| quality    | Evaluate the quality of your Hi-C data                                                                            |
+------------+-------------------------------------------------------------------------------------------------------------------+

Links
=====
- `Detailed Documentation <http://xiaotaowang.github.io/HiC_pipeline/>`_
    - `Installation <http://xiaotaowang.github.io/HiC_pipeline/install.html>`_
    - `Quick Start <http://xiaotaowang.github.io/HiC_pipeline/quickstart.html>`_
    - `Data Quality <http://xiaotaowang.github.io/HiC_pipeline/quality.html>`_
    - `Parallel Computation <http://xiaotaowang.github.io/HiC_pipeline/parallel.html>`_
- `Code Repository <https://github.com/XiaoTaoWang/HiC_pipeline/>`_ (At GitHub, Track the package issue)
- `PyPI <https://pypi.python.org/pypi/runHiC>`_ (Download and Installation)

Usage
=====
Open a terminal, type ``runHiC -h`` or ``runHiC <subcommand> -h`` for help information.

Citation
========
Xiaotao Wang. (2016). runHiC: A user-friendly Hi-C data processing software based on hiclib. Zenodo.
`10.5281/zenodo.55324 <http://dx.doi.org/10.5281/zenodo.55324>`_
