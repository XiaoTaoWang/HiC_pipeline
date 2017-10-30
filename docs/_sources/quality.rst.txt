Experiment Quality Assessment
=============================
In this tutorial, I will show you how runHiC can be used in data quality assessment.

All you need to type in is a sinlge line command after *runHiC filtering* or *runHiC pileup*
(refer to `quickstart <http://xiaotaowang.github.io/HiC_pipeline/quickstart.html>`_ for more details)::

    $ runHiC quality -m datasets.tsv -L filtered-hg19

- ``-m/--metadata``

  The metadata data file name.
	
- ``-L/--Locator``

  Path to the folder containing filtered HDF5 files.

Statistic Table
---------------
In our example, statistic tables on sequencing reads for each SRA/FASTQ (level 1),
biological replicate (level 2) and cell type (level 3) will be generated under the
"filtered-hg19" folder.

Here's a snapshot::

    000_SequencedReads:   14332993
	    010_UniqueMappedReads:   13053095
	    020_LigationCounts:   3225524
    100_DoubleUniqueMapped:   8260823
	    110_AfterFilteringReads:   6836510
	    120_SameFragmentReads:   1147363
		    122_SelfLigationReads:   210202
		    124_DanglingReads:   928888
		    126_UnknownMechanism:   8273
	    210_ExtraDanglingReads:   250193
	    310_DuplicatedRemoved:   26757
    400_TotalContacts:   6836510
	    410_IntraChromosomalReads:   2753015
		    412_IntraLongRangeReads(>=20Kb):   2402139
		    412_IntraShortRangeReads(<20Kb):   350876
	    420_InterChromosomalReads:   4083495

    Critical Indicators:
    Double Unique Mapped Ratio = 8260823 / 14332993 = 0.5764
    Ligation-Junction Ratio = 3225524 / 14332993 = 0.2250
    Self-Ligation Ratio = 210202 / 14332993 = 0.0147
    Dangling-Reads Ratio = 928888 / 14332993 = 0.0648
    Long-Range Ratio = 2402139 / 6836510 = 0.3514
    Data Usage = 6836510 / 14332993 = 0.4770

The following table lists possible statistic names and their meanings:

+-------------------------------+---------------------------------------------------+
| Statistic Name                | Meaning                                           |
+===============================+===================================================+
| 000_SequencedReads            | Total number of sequenced read pairs              |
+-------------------------------+---------------------------------------------------+
| 010_UniqueMappedReads         | Number of read pairs of which one read or both    |
|                               | reads can be uniquely mapped to the reference     |
|                               | genome.                                           |
+-------------------------------+---------------------------------------------------+
| 020_LigationCounts            | Number of read pairs containing the so-called     |
|                               | "ligation junction". A ligation junction is       |
|                               | the sequence created when the ends of two         |
|                               | filled-in restriction fragments ligate to one     |
|			        | another. For HindIII, the sequence is AAGCTAGCTT, |
|                               | and for MboI, it's GATCGATC. Obviously, this      |
|                               | statistic is dependent on the read length and     |
|                               | your library size. For 300~500bp library and      |
|                               | 101bp PE reads, it generally falls into the       |
|                               | 30%~40% range. A low value suggests that the      |
|                               | ligation failed. [1]_                             |
+-------------------------------+---------------------------------------------------+
| 100_DoubleUniqueMapped        | Number of read pairs of which both reads can be   |
|                               | uniquely mapped.                                  |
+-------------------------------+---------------------------------------------------+
| 110_AfterFilteringReads       | Number of read pairs that have passed all         |
|                               | filtering criteria.                               |
+-------------------------------+---------------------------------------------------+
| 120_SameFragmentReads         | Number of read pairs of which both reads are      |
|                               | mapped to the same restriction fragment. Such     |
|                               | read pairs are filtered in our pipeline.          |
+-------------------------------+---------------------------------------------------+
| 122_SelfLigationReads         | Number of read pairs deriving from                |
|                               | self-circularized ligation product. The two reads |
|                               | are mapped to the same restriction fragment and   |
|                               | face in opposite directions.                      |
+-------------------------------+---------------------------------------------------+
| 124_DanglingReads             | Both reads of these read pairs are mapped to the  |
|                               | same fragment and face toward each other. There   |
|                               | can be many causes of such products, ranging from |
|                               | low ligation efficiency to poor streptavidin      |
|                               | specificity.                                      |
+-------------------------------+---------------------------------------------------+
| 126_UnknownMechanism          | Unknown sources of "120_SameFragmentReads". Both  |
|                               | reads are mapped to the same strand.              |
+-------------------------------+---------------------------------------------------+
| 210_ExtraDanglingReads        | The two reads of these read pairs are mapped to   |
|                               | different restriction fragments but face toward   |
|                               | each other and are separated by less than the     |
|                               | library size (500bp) interval. Such read pairs    |
|                               | may contain true contacts, but are largely        |
|                               | contaminated, so we also remove these read pairs  |
|                               | from our analysis.                                |
+-------------------------------+---------------------------------------------------+
| 310_DuplicatedRemoved         | Number of read pairs from PCR products. We treat  |
|                               | two read pairs to be duplicated from one another  |
|                               | if both reads of them are mapped to the same      |
|                               | position of the genome. Such redundant read pairs |
|                               | are also filtered from our analysis.              |
+-------------------------------+---------------------------------------------------+
| 400_TotalContacts             | Number of read pairs from true contacts, i.e.,    |
|                               | the remaining read pairs after all filtering      |
|                               | processes.                                        |
+-------------------------------+---------------------------------------------------+
| 410_IntraChromosomalReads     | Number of intra-chromosomal contacts              |
+-------------------------------+---------------------------------------------------+
| 412_IntraLongRangeReads       | Number of long-range contacts (genomic distance   |
|                               | >= 20Kb)                                          |
+-------------------------------+---------------------------------------------------+
| 412_IntraShortRangeReads      | Number of short-range contacts (genomic distance  |
|                               | < 20Kb)                                           |
+-------------------------------+---------------------------------------------------+
| 420_InterChromosomalReads     | Number of inter-chromosomal contacts              |
+-------------------------------+---------------------------------------------------+
| 500_IntraMitochondrial        | Number of intra-mitochondrial contacts            |
+-------------------------------+---------------------------------------------------+
| 600_InterNuclearMitochondrial | Number of contacts between mitochondrial genome   |
|                               | and the nuclear genome. This indicator has        |
|                               | potential to assess the random ligation level of  |
|                               | your library.                                     |
+-------------------------------+---------------------------------------------------+

Note that we try to organize these statistics hierarchically using indentation,
so "010_UniqueMappedReads" and "020_LigationCounts" are parts of "000_SequencedReads",
similarly, "122_SelfLigationReads", "124_DanglingReads" and "126_UnknownMechanism"
constitute "120_SameFragmentReads".

At the bottom of the statistic table, we include some important quality indicators:

1. Unique-Mapping Ratio. Low value of this metric indicates low sequencing quality,
   sample contamination or incomplete genome assembly.

2. Ligation-Junction Ratio. Low value of this metric suggests the failure of ligation.

3. Self-Ligation Ratio.

4. Dangling-Reads Ratio.

5. Long-Range Ratio. Low value (<0.15) of this metric indicates the failed experiment.

Library-size Estimation
------------------------
Dangling reads can be applied to estimate your library size in nature. Here's an example
of size distribution of dangling read molecules for typical 300~500bp library:

.. image:: ./_static/GM06990-HindIII-allReps-librarySize.png
        :align: center

The inconsistency between this distribution and the experimental library size suggests
a failure in DNA size selection.

Ligation Efficiency
-------------------
Excessive dangling reads may result from low ligation efficiency or poor streptavidin
specificity. Can we further discriminate these two mechanisms? Intuitively, if one of
reads of most read pairs locate near a restriction site, the former (low ligation efficiency)
is more likely to be the cause, so we also plot the distribution of the relative start
sites for dangling reads:

.. image:: ./_static/GM06990-HindIII-allReps-danglingStart.png
        :align: center

Here, the majority of these read pairs have one of their read starting near a restriction
site, therefore, ligation efficiency could be a good explain.

Read-pair Type Plotting
-----------------------
Read-pair type ratios will be reported in line-plot manner for each biological
replicate (level 1) and cell type (level 2) under filtered-hg19 too. Intra-chromosomal
contacts are broken down into four types: "left pair" (both reads map to the reverse
strand), "right pair" (both reads map to the forward strand), "inner pair" (reads map
to different strands and point towards each other) and "outer pair" (reads map to
different strands and point away from one another). If reads come from proximity
ligation, each pair type should account for roughly 25% of contacts. Thus, distance
at which the percentage of each type converges to 25% is a good indication of the minimum
distance at which it is meaningful to examine Hi-C contact patterns. Here's an example
below:

.. image:: ./_static/GM06990-HindIII-allReps-PairType.png
        :align: center

We can see a distinct turning point around 20Kb. While there may be several unknown mechanisms
making biases below this point, we should only consider contacts whose genomic distances
are greater than 20Kb in the following analysis.


Reference
=========
.. [1] Rao SS, Huntley MH, Durand NC et al. A 3D Map of the Human Genome at Kilobase Resolution
       Reveals Principles of Chromatin Looping. Cell, 2014, 159(7):1665-80.
