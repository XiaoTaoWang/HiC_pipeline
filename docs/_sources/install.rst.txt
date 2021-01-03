Install runHiC
==============

Requirements
------------
runHiC is developed and tested on UNIX-like operating system, and following packages
or softwares are required:

Python requirements:

- Python 3.7+
- numpy 1.19+
- matplotlib 2.0+
- biopython 1.78+
- pairtools 0.3.0+
- cooler 0.8.6

Other requirements:

- sra-tools 2.10.8+
- bwa 0.7.17+
- minimap2 2.17+
- samtools 1.10+

Optional:

- pigz 2.3.4+

Install Requirements through Conda
----------------------------------
All above requirements can be installed through the conda package manager.

.. note:: If you have the Anaconda Distribution installed, you already have it.

Choose an appropriate `Miniconda installer <https://conda.io/miniconda.html>`_ for your system,
then in your terminal window, type the following and follow the prompts on the installer screens::

    $ bash Miniconda3-latest-Linux-x86_64.sh

After that, update the environment variables to finish the Conda installation::

    $ source ~/.bashrc

Conda allows separation of packages into separate repositories, or channels. The main *defaults*
channel only covers *numpy*, *matplotlib* and *biopython* listed above. And to make all the other packages
accessible, you need to add the *bioconda* channel and the *conda-forge* channel in the following way (note
that the order is important to guarantee the correct priority)::

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge

Then type and execute the commands below to satisfy the requirements::

    conda create -n runHiC numpy matplotlib biopython pairtools cooler=0.8.6 sra-tools bwa minimap2 samtools pigz
    conda activate runHiC

Install runHiC
--------------
Finally, *runHiC* can be installed from PyPI using pip::
    
    $ pip install runHiC

runHiC has been installed successfully if no exception occurs in the above processes.
