Install runHiC
==============

Requirements
------------
runHiC is developed and tested on UNIX-like operating system, and following packages
or softwares are required:

Python requirements:

- Python (2.7, not compatible with 3.x for now)
- numpy
- numexpr 
- scipy
- statsmodels
- matplotlib
- h5py
- cython
- biopython
- bx-python
- pysam
- joblib
- pp
- mirnylib
- hiclib

Other requirements:

- gcc (>= 4.8.5)
- sra-tools
- bowtie2
- samtools

Optional:

- pigz

Install Conda
-------------
All above requirements except for mirnylib and hiclib can be installed through the
conda package manager.

.. note:: If you have the Anaconda Distribution installed, you already have it, feel free to jump to
   the `Set up Channels`_ section.

Download the latest `Linux Miniconda installer for Python 2.7 <https://conda.io/miniconda.html>`_,
then in your terminal window type the following and follow the prompts on the installer screens::

    $ bash Miniconda2-latest-Linux-x86_64.sh

After that, update the environment variables to finish the Conda installation::

    $ source ~/.bashrc

Set up Channels
---------------
Conda allows separation of packages into separate repositories, or channels. The main *defaults*
channel has a large amount of common packages including *numpy*, *numexpr*, *scipy*, *statsmodels*,
*matplotlib*, *h5py*, *cython* and *biopython* listed above. *bx-python*, *pysam*, *joblib*, *pp*,
*sra-tools*, *bowtie2* and *samtools* are not available in the *defaults* channel but included in
the *bioconda* channel, and to make them accessible, you will need to add the *bioconda* channel
as well as the other channels bioconda depends on (note that the order is important to guarantee
the correct priority)::

    $ conda config --add channels conda-forge
    $ conda config --add channels defaults
    $ conda config --add channels r
    $ conda config --add channels bioconda

Install Packages through Conda
------------------------------
Then it's straightforward to install all the required packages except for *mirnylib* and *hiclib*
through the following one-line command::

    conda install numpy numexpr scipy statsmodels matplotlib h5py cython biopython bx-python pysam joblib pp sra-tools bowtie2 samtools pigz gcc=4.8.5

Install mirnylib and hiclib
---------------------------
To make it easy, I have included the *mirnylib* and *hiclib* source code in the `runHiC <https://pypi.python.org/pypi/runHiC>`_
distribution under the "mirnylab" sub-foler since version 0.7.0::

    $ cd mirnylab

To install *mirnylib*::

    $ unzip mirnylib.zip
    $ cd mirnylib
    $ python setup.py install

To install *hiclib*::

    $ unzip hiclib.zip
    $ cd hiclib
    $ python setup.py install

Install runHiC
--------------
Now just run the setup.py script under the distribution root folder to finish the installation.
Suppose you are still in the extracted hiclib folder::

    $ cd ../..
    $ python setup.py install

runHiC has been installed successfully if no exception occurs in the above processes.
