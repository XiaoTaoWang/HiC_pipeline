runHiC
******
.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.55324.svg
   :target: http://dx.doi.org/10.5281/zenodo.55324

runHiC is an easy-to-use command-line tool for Hi-C data processing.

Since version 0.8.1, runHiC can be used directly on `Arima HiC <https://arimagenomics.com>`_ data
by setting the enzyme name to *Arima*.

Since version 0.8.0, runHiC has changed its default data container/format from HDF5 to
`Pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_ and
`Cooler <https://github.com/mirnylab/cooler>`_. (See `Release Notes <http://xiaotaowang.github.io/HiC_pipeline/changelog.html>`_)

Design Concepts
===============
runHiC is designed to process Hi-C data from raw sequencing reads(.sra, .fastq, .fastq.gz) to the corrected
contact matrices. It currently contains 5 subcommand:

+------------+-------------------------------------------------------------------------------------+
| mapping    | Map raw pair-end sequencing data to a supplied genome. Support bwa and minimap2.    |
+------------+-------------------------------------------------------------------------------------+
| filtering  | Perform read-level and fragment-level noise removing                                |
+------------+-------------------------------------------------------------------------------------+
| binning    | 1.Generate contact matirx; 2. Perform ICE                                           |
+------------+-------------------------------------------------------------------------------------+
| pileup     | Perform entire processing from *mapping* to *binning*                               |
+------------+-------------------------------------------------------------------------------------+
| quality    | Assess the quality of your Hi-C data                                                |
+------------+-------------------------------------------------------------------------------------+

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
Open a terminal, type ``runHiC -h`` and ``runHiC <subcommand> -h`` for help information.

Citation
========
Xiaotao Wang. (2016). runHiC: A user-friendly Hi-C data processing software based on hiclib. Zenodo.
`10.5281/zenodo.55324 <http://dx.doi.org/10.5281/zenodo.55324>`_
