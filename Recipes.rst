Suppose your directory structure looks exactly same as the **Sample** folder
distributed with our source code.

Basic Pipeline
--------------
Change to the working directory::

    $ cd Sample/working

Simply run following command in a terminal prompt::

    $ runHiC pileup -p ../data -g hg19 --fastqDir SRA -F SRA --bowtiePath ../Tools/bowtie2/bowtie2 -m datasets.tsv --chunkSize 1500000 --libSize 500

*pileup* streamlines all analysis stages from mapping to ICE correcting.

If you want to perform a step-by-step analysis, you need to call ``mapping``,
``filtering``, ``binning``, ``correcting`` and ``tosparse`` in order.

Mapping::

    $ runHiC mapping -p ../data -g hg19 --fastqDir SRA -F SRA --bowtiePath ../Tools/bowtie2/bowtie2 -m datasets.tsv --chunkSize 1500000

Filtering::

    $ runHiC filtering -p ../data -g hg19 --HDF5 hdf5-hg19 -m datasets.tsv --libSize 500 --duplicates --startNearRsite --level 2

Binning::

    $ runHiC binning -p ../data -g hg19 --filteredDir filtered-hg19 --mode wholeGenome --resolution 200000

Correcting::

    $ runHiC correcting -p ../data -g hg19 --HeatMap Heatmaps-hg19
	
Convert to sparse format::

    $ runHiC tosparse -p ../data -g hg19 --cHeatMap Corrected-hg19

