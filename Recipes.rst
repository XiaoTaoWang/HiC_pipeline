Suppose your directory structure looks exactly same as the **Sample** folder
distributed with our source code.

Basic Pipeline
--------------
Change to the working directory::

    $ cd Sample/working

Simply run following command in a terminal prompt::

    $ runHiC pileup -p ../data -g hg19 --fastqDir SRA -F SRA --bowtieIndex ../Tools/bowtie2/bowtie2 -m datasets.tsv --chunkSize 1500000

*pileup* streamlines all analysis stages from mapping to ICE correcting.