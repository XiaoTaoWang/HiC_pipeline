Suppose your directory structure looks exactly same as the **Sample** folder
distributed with our source code.

Basic Pipeline
--------------
Change to the working directory::

    $ cd Sample/working

Simply run following command in a terminal prompt::

    $ runHiC pileup -p ../data -g hg19 --fastqDir SRA -F SRA --bowtiePath ../Tools/bowtie2/bowtie2 -m datasets.tsv --chunkSize 1500000 --libSize 500

*pileup* streamlines all analysis stages from mapping to ICE correcting.

If you want to perform a step-by-step analysis, you need to call *mapping*,
*filtering*, *binning*, *correcting* and *tosparse* separately, in order.

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

Parallel Tasks
--------------
a) Bowtie2 supports multiple threads for alignments. (You can specify the number
   of threads through ``-t/--threads`` when running *mapping* or *pileup*)
b) runHiC provides another layer for parallel computing. On this level of parallel,
   tasks are arranged based on separate SRA/FASTQ file, i.e., you can use this
   capacity only if you have two or more SRA/FASTQ files. Just submit the same command
   repeatedly and the program allocate a unique ID for each SRA/FASTQ to avoid conflicts
   between processes.

I give an example below to use these capacities as much as possible.

At first, run the similar but slightly changed command below::

    $ runHiC pileup -p ../data -g hg19 --fastqDir SRA -F SRA --bowtiePath ../Tools/bowtie2/bowtie2 -t 10 --logFile runHiC-pileup-1.log -m datasets.tsv --chunkSize 1500000 --libSize 500
	
This command uses 10 bowtie2 threads and redirect logging messages to another file
named "runHiC-pileup-1.log".

Several seconds later, run (If you are computing on a cluster, you may need to
switch to another node for efficiency)::

    $ runHiC pileup -p ../data -g hg19 --fastqDir SRA -F SRA --bowtiePath ../Tools/bowtie2/bowtie2 -t 10 --logFile runHiC-pileup-2.log -m datasets.tsv --chunkSize 1500000 --libSize 500
	
This time, logging messages are written to "runHiC-pileup-2.log".

That's not all, if you have three or more SRA files::

    $ runHiC pileup -p ../data -g hg19 --fastqDir SRA -F SRA --bowtiePath ../Tools/bowtie2/bowtie2 -t 10 --logFile runHiC-pileup-3.log -m datasets.tsv --chunkSize 1500000 --libSize 500
	
...

Experiment Quality Assessment
-----------------------------
Call *quality* after *filtering* or *pileup*::

    $ runHiC quality -p ../data -g hg19 -m datasets.tsv

Statistic table on sequencing reads for each SRA/FASTQ (level 1), biological
replicate (level 2) and cell type (level 3) will be generated under filtered-hg19.

Read-pair type ratios will be reported in line-plot manner for each biological
replicate (level 2) and cell type (level 3) under filtered-hg19 too. Intra-chromosomal
contacts are broken down into four types: "left" pair (both reads map to the reverse
strand), "right" pair (both reads map to the forward strand), "inner" pair (reads map
to different strands and point towards each other) and "outer" pair (reads map to
different strands and point away from one another). If the reads come from proximity
ligation, each pair type should account for roughly 25% of contacts. Thus, distance
at which the percentage of each type converges to 25% is a good indication of the minimum
distance at which it is meaningful to examine Hi-C contact patterns.
