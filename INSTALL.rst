Python Packages
===============
I recommend using `conda <http://conda.pydata.org/miniconda.html>`_, an excellent Python package and
environment manager.

Open a terminal and type::

    $ conda install numpy numexpr scipy matplotlib cython biopython h5py pysam pip

Install *bx-python* and *joblib* using pip::

    $ pip install joblib bx-python

Install *mirnylib* and *hiclib* from source code:

Download `mirnylib <https://bitbucket.org/mirnylab/mirnylib>`_ and `hiclib <https://bitbucket.org/mirnylab/hiclib>`_,
and run install_linux.py contained in the unpacked folder, respectively.

Other dependencies
==================
Install *fastq-dump*:

I have included it in our distribution. Just make it accessible to your system. (Via the environment variable
*PATH*)

Install *samtools*:

Download `samtools <http://sourceforge.net/projects/samtools/files/>`_, unpack it, change to the extracted
directory::

    $ make

Make *samtools* accessible to your system.

Install *Bowtie2*:

Download the `source code <http://sourceforge.net/projects/bowtie-bio/files/bowtie2/>`_, unzip it and
add the path to the extracted directory to *PATH*.

runHiC
======
Use easy_install::

    $ easy_install runHiC