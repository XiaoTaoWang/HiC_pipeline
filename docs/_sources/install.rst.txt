Install runHiC
==============

Requirements
------------
runHiC is developed and tested on UNIX-like operating systems, and following packages
or softwares are required:

Python requirements:

- Python 3.7+
- matplotlib 2.0+
- biopython 1.78+
- pairtools 0.3.0+
- cooler 0.8.6+

Other requirements:

- sra-tools 2.10.8+
- bwa 0.7.17+
- minimap2 2.17+
- chromap 0.1+
- samtools 1.10+
- pigz 2.3.4+

Install Requirements through Mamba
----------------------------------
To satisfy the requirements above, just install `mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_
on your machine, and execute the commands below::

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ mamba create -n runHiC matplotlib biopython pairtools cooler sra-tools bwa minimap2 samtools pigz chromap
    $ mamba activate runHiC

Install runHiC
--------------
Then *runHiC* can be installed from PyPI using pip::
    
    $ pip install runHiC

You should have runHiC installed successfully if no exception occurs in the above processes.
