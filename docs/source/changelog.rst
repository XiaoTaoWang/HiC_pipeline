Release Notes
=============
Version 0.8.7 (01/19/2023)
--------------------------
- Made the program compatible with pairtools v1.0.2

Version 0.8.6 (05/24/2022)
--------------------------
- Added support for other 3C-based experiments that do no use restriction enzymes for DNA fragmentation

Version 0.8.5 (08/15/2021)
---------------------------
- changed the default read aligner to chromap
- made the fragment-level filtering optional
- removed the "-O/--outformat" argument

Version 0.8.4-r1 (01/16/2021)
------------------------------
- Added the "--high-res" option to ``runHiC pileup`` and ``runHiC binning``.

Version 0.8.4 (01/03/2021)
--------------------------
- Changed the bwa parameters from "bwa mem -SP" to "bwa mem -SP5M"
- Changed the default walks policy from "mask" to "all"
- Fixed a bug pertaining to the FASTQ read format. In previous versions, error will occur when there is "." in the prefix of FASTQ file names.
- Added a piechart function to ``runHiC quality``.
- Added a function to automatically output mcool in the final step.
- Added a function to automatically check and remove invalid link files within the "filtered" sub-folder
- Fixed a bug to avoid duplicate operations when only one replicate exists
- Added a function to automatically sort chromosome names. Numerical labels will be sorted natually; for non-numerical labels, give priority to XYM over the rest.

Version 0.8.3 (03/04/2019)
--------------------------
- Parallelized statistics collecting
- Reduced the computational overhead in the filtering step

Version 0.8.2 (01/16/2019)
--------------------------
- Fixed a bug pertaining to large dataset overflow in binning (by increasing the value for ``--max-split``)

Version 0.8.1 (12/23/2018)
--------------------------
- Added support for Arima Hi-C

Version 0.8.0 (09/16/2018)
--------------------------
- Migrated to Python 3 (will not support Python 2 anymore)
- Changed the data container from HDF5 to Pairs and Cooler
- Added two aligner options: bwa-mem (short reads) / minimap2 (long reads)
- Added options to output filtered SAM/BAM for further analysis

Version 0.7.0 (10/29/2017)
---------------------------
- Parallelized the execution on both laptops and PBS-based clusters.
- Included the latest *mirnylib* and *hiclib* source code.
- Tuned the class/function interfaces in sync with the lastest *hiclib*.
- Removed the *visualize* subcommand.
- Changed the default "--chroms" setting to "['#','X']"
- Refined the logging information.
- Added the detailed documentations.

Version 0.6.6-r2 (05/04/2016)
-----------------------------
- Fixed a bug related to "--bowtieIndex".
- Included LICENSE.

Version 0.6.6-r1 (11/10/2015)
-----------------------------
- Only check for update when you are online.

Version 0.6.6 (11/06/2015)
--------------------------
- Removed redundant command-line arguments.
- Added detailed descriptions for the quality assessment module in Recipes.rst.
- Added snippets for automatic update checking.

Version 0.6.5 (09/29/2015)
--------------------------
- More sophisticated quality assessment module.
- Changed the filtering settings for *pileup*.
- Redirected the exception information into the log file.

Version 0.6.4 (09/21/2015)
--------------------------
- Fixed a bug related to the metadata loading.

Version 0.6.3 (09/16/2015)
--------------------------
- Customized the mirnylib Genome class by overriding the *_extractChrmLabel* method.

Version 0.6.2 (08/21/2015)
--------------------------
- Added the *visualize* subcommand.

Version 0.6.1 (07/21/2015)
--------------------------
- Improved the logging system.

Version 0.6.0 (06/16/2015)
--------------------------
- Fixed a bug for bam file parsing in the case of FASTQ read format.
