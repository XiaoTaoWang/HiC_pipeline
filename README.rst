runHiC
******
.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.55324.svg
   :target: http://dx.doi.org/10.5281/zenodo.55324

runHiC is a easy-to-use Hi-C processing software based on hiclib (https://bitbucket.org/mirnylab/hiclib)
Different from hiclib, which was born for flexibility, runHiC is a customized pipeline, and can be
run from command line directly.

Since version 0.7.0, runHiC has been able to execute concurrently either on a single laptop with multiple
processors or on a PBS-based clusters.

Design Concepts
===============
runHiC is designed to process Hi-C data from raw sequencing reads(.sra, .fastq, .fastq.gz) to the corrected
contact matrices. It separates the whole procedure into 4 stages(*mapping*, *filtering*, *binning*,
*correcting*) and contains 7 subcommands:

+------------+------------------------------------------------------------------------------+
| mapping    | Iteratively map pair-end sequencing reads to a supplied genome               |
+------------+------------------------------------------------------------------------------+
| filtering  | Perform read-level and fragment-level noise removing                         |
+------------+------------------------------------------------------------------------------+
| binning    | Generate the original contact matrices                                       |
+------------+------------------------------------------------------------------------------+
| correcting | Perform iterative corrections on original contact matrices                   |
+------------+------------------------------------------------------------------------------+
| tosparse   | Convert the dense intra-chromosomal contact matrices to sparse ones          |
+------------+------------------------------------------------------------------------------+
| pileup     | Streamline all stages from *mapping* to *correcting*                         |
+------------+------------------------------------------------------------------------------+
| quality    | Assess the quality of your Hi-C data                                         |
+------------+------------------------------------------------------------------------------+

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
