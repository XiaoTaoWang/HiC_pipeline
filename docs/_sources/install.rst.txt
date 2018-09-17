Install runHiC
==============

Requirements
------------
runHiC is developed and tested on UNIX-like operating system, and following packages
or softwares are required:

Python requirements:

- Python 3.5+
- numpy
- matplotlib
- biopython
- pairtools
- cooler

Other requirements:

- sra-tools
- bwa
- minimap2
- samtools

Optional:

- pigz

Install Requirements through Conda
----------------------------------
All above requirements can be installed through the conda package manager.

.. note:: If you have the Anaconda Distribution installed, you already have it.

Choose an appropriate `Miniconda installer <https://conda.io/miniconda.html>`_ for your system,
then in your terminal window type the following and follow the prompts on the installer screens::

    $ bash Miniconda3-latest-Linux-x86_64.sh

After that, update the environment variables to finish the Conda installation::

    $ source ~/.bashrc

Conda allows separation of packages into separate repositories, or channels. The main *defaults*
channel only covers *numpy*, *matplotlib* and *biopython* listed above. And to make all the other packages
accessible, you need to add the *bioconda* channel and *conda-forge* channel in the following way (note
that the order is important to guarantee the correct priority)::

    $ conda config --add channels conda-forge
    $ conda config --add channels defaults
    $ conda config --add channels r
    $ conda config --add channels bioconda

Then just type and execute this one-line command to satisfy the requirements::

    conda install numpy matplotlib biopython pairtools cooler sra-tools bwa minimap2 samtools pigz

Install runHiC
--------------
Finally, *runHiC* can be installed from PyPI by pip::
    
    $ pip install runHiC

runHiC has been installed successfully if no exception occurs in the above processes.
