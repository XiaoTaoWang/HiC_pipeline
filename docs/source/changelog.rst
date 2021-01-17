Release Notes
=============
Version 0.8.4-r1 (2021-01-16)
-----------------------------
- added the "--high-res" option to ``runHiC pileup`` and ``runHiC binning``.

Version 0.8.4 (2021-01-03)
--------------------------
- Changed the bwa parameters from "bwa mem -SP" to "bwa mem -SP5M"
- Updated the default walks policy from "mask" to "all"
- Added the option "--memory" to change the allocated memory for ``pairtools sort`` and ``pairtools merge``.
- Fixed a bug pertaining to the FASTQ read format. Previously, error occurred when there is "." in the prefix of FASTQ file names.
- Added a piechart function to plot the percentage of each category for ``runHiC quality``.
- Automatically output mcool (1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000) in the final step.
- Automatically check and remove invalid link files within the "filtered" folder
- Avoid duplicate operations when only one replicate exists
- Automatically sort chromosome names. Numerical labels will be sorted natually; for non-numerical labels, give priority to XYM over the rest.

Version 0.8.3 (2019-03-04)
--------------------------
- Parallelized statistics collecting
- Reduced the computational overhead in the filtering step

Version 0.8.2 (2019-01-16)
--------------------------
- Fixed the bug of large dataset overflow in binning (by increasing the value for ``--max-split``)

Version 0.8.1 (2018-12-23)
--------------------------
- Supported Arima Hi-C

Version 0.8.0 (2018-9-16)
-------------------------
- Migrated to Python 3 (will not support Python 2 anymore)
- Changed the data container from HDF5 to Pairs and Cooler
- Added two aligner options: bwa-mem (short reads) / minimap2 (long reads)
- Added options to output filtered SAM/BAM for further analysis

Version 0.7.0 (2017-10-29)
--------------------------
- Parallelized the execution on both laptops and PBS-based clusters.
- Included the latest *mirnylib* and *hiclib* source code.
- Tuned the class/function interfaces in sync with the lastest *hiclib*.
- Removed the *visualize* subcommand.
- Changed the default "--chroms" setting to "['#','X']"
- Refined the logging information.
- Added the detailed documentations.

Version 0.6.6-r2 (2016-5-4)
---------------------------
- Fixed a bug related to "--bowtieIndex".
- Included LICENSE.

Version 0.6.6-r1 (2015-11-10)
-----------------------------
- Only check for update when you are online.

Version 0.6.6 (2015-11-6)
-------------------------
- Removed redundant command-line arguments.
- Added detailed descriptions for the quality assessment module in Recipes.rst.
- Added snippets for automatic update checking.

Version 0.6.5 (2015-9-29)
-------------------------
- More sophisticated quality assessment module.
- Changed the filtering settings for *pileup*.
- Redirected the exception information into the log file.

Version 0.6.4 (2015-9-21)
-------------------------
- Fixed a bug related to the metadata loading.

Version 0.6.3 (2015-9-16)
-------------------------
- Customized the mirnylib Genome class by overriding the *_extractChrmLabel* method.

Version 0.6.2 (2015-8-21)
-------------------------
- Added the *visualize* subcommand.

Version 0.6.1 (2015-7-21)
-------------------------
- Improved the logging system.

Version 0.6.0 (2015-6-16)
-------------------------
- Fixed a bug for bam file parsing in the case of FASTQ read format.
