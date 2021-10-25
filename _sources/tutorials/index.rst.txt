.. _tutorials:

Tutorials
=========

^^^^^

BioKIT can be used for a multitude of different types of analyses. Documentation here 
provides a step-by-step outline for how to conduct some of these analyses.

1. summarize genome assembly metrics,
2. summarize multiple sequence alignment properties,
3. calculate relative synonymous codon usage (RSCU), and 
4. calculate a novel metric for estimating codon optimization, gene-wise RSCU (or gw-RSCU). 

|

1. Summarize genome assembly metrics
####################################

^^^^^

BioKIT implements numerous functions that can be used to evaluate diverse properties of a genome assembly.

In this tutorial, we will examine the contig-level genome assembly of *Escherichia coli* strain ATCC 11775,
which was obtained from NCBI in September of 2021. |br|

.. centered::
   Download test data:
   :download:`E. coli strain ATCC 11775 genome assembly</data/GCA_000734955.1_ASM73495v1_genomic.fna>`

To calculate all genome assembly metrics, use the following command:

.. code-block:: shell

   biokit genome_assembly_metrics GCA_000734955.1_ASM73495v1_genomic.fna 
   
   4980585 Assembly size
   14      L50
   40      L90
   124187  N50
   33720   N90
   0.5062  GC content
   135     Number of scaffolds
   93      Number of large scaffolds
   4968723 Sum length of large scaffolds
   382114  Longest scaffold
   0.247   Frequency of A
   0.2468  Frequency of T
   0.2525  Frequency of C
   0.2537  Frequency of G

If all metrics are not needed, individual metrics can be calculated. For example,
calculate N50 of the genome assembly using the following command:

.. code-block:: shell

   biokit n50 GCA_000734955.1_ASM73495v1_genomic.fna 
   
   124187

|

2. Summarize multiple sequence alignment properties
###################################################

^^^^^

The properties of multiple sequence alignments are commonly reported in studies. To
demonstrate how to calculate diverse properties of multiple sequence alignments, we
will use the following single gene alignment from `Steenwyk et al. 2019, mBio
<https://journals.asm.org/doi/10.1128/mbio.00925-19?permanently=true>`_. 

.. centered::
   Download test data:
   :download:`Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa</data/Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa>`


To calculate diverse properties of the multiple sequence alignment, use the following
command:

.. code-block:: shell

   biokit alignment_summary Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa 
   
   General Characteristics
   =======================
   90      Number of taxa
   624     Alignment length
   517     Parsimony informative sites
   555     Variable sites
   69      Constant sites

   Character Frequencies
   =====================
   T       8527
   G       14835
   C       15482
   A       12966
   -       4350

If all metrics are not needed, individual metrics can be calculated. For example,
the number of constant sites can be calculated using the following command:

.. code-block:: shell

   biokit constant_sites Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa 
   
   69

|

3. Calculate relative synonymous codon usage bias
#################################################

^^^^^

Patterns of codon usage can be examined to reveal biases favoring some codons over others.
In this section, we will calculate patterns of relative synonymous codon usage (RSCU) among
genes encoded in the genome of *Saccharomyces cerevisiae*, the model budding yeast.

.. centered::
   Download test data:
   :download:`GCA_000146045.2_R64_cds_from_genomic.fna</data/GCA_000146045.2_R64_cds_from_genomic.fna>`

To calculate RSCU, use the following command:

.. code-block:: shell

   biokit rscu GCA_000146045.2_R64_cds_from_genomic.fna
   
   AGA     2.8425
   GGU     1.8179
   UUG     1.6747
   UUA     1.6645
   CCA     1.6281
   UCU     1.5608
   GUU     1.5466
   GCU     1.4773
   UAA     1.4221
   GAA     1.4014
   AUU     1.3807
   CAA     1.371
   ACU     1.3705
   GAU     1.3024
   CAU     1.2838
   AGG     1.2789
   UCA     1.2736
   UGU     1.2424
   CCU     1.2412
   ACA     1.2315
   AAU     1.1923
   UUU     1.1892
   GCA     1.1865
   AAA     1.1684
   UAU     1.1331
   AUG     1.0
   UGG     1.0
   AGU     0.9759
   UCC     0.9403
   GGA     0.903
   UGA     0.888
   GCC     0.8851
   GUA     0.8705
   UAC     0.8669
   CUA     0.8549
   CGU     0.849
   ACC     0.8452
   AUA     0.8406
   AAG     0.8316
   UUC     0.8108
   GUC     0.8093
   AAC     0.8077
   GGC     0.7898
   CUU     0.7804
   AUC     0.7787
   GUG     0.7736
   UGC     0.7576
   CAC     0.7162
   GAC     0.6976
   UAG     0.6899
   CUG     0.6745
   AGC     0.6667
   CCC     0.6332
   CAG     0.629
   GAG     0.5986
   UCG     0.5827
   ACG     0.5528
   CCG     0.4974
   GGG     0.4892
   GCG     0.4512
   CGA     0.4217
   CGC     0.3587
   CUC     0.351
   CGG     0.2491

Values greater than 1 are reflective of codons that are over represented whereas
values less than 1 are reflective of codons that are under represented and values 
equal to one occur as expected under a neutral expectation. Here, AGA is the most
over represented codon and CGG is the most under represented codon.

|

4. Estimating codon optimization
################################

^^^^^

Estimates of codon optimization can be used to identify genes that are expressed at 
high levels and are important for organismal ecology. Here, we will calculate gw-RSCU,
a measure of codon optimization, in the coding sequences of *Saccharomyces cerevisiae*,
the model budding yeast.

.. centered::
   Download test data:
   :download:`GCA_000146045.2_R64_cds_from_genomic.fna</data/GCA_000146045.2_R64_cds_from_genomic.fna>`

To calculate gw-RSCU, use the following command:

.. code-block:: shell

   biokit gw_rscu GCA_000146045.2_R64_cds_from_genomic.fna 

   lcl|BK006935.2_cds_DAA06918.1_1 1.1972  1.0     0.4678
   lcl|BK006935.2_cds_DAA06919.1_2 1.0037  0.9515  0.3458
   lcl|BK006935.2_cds_DAA06920.1_3 1.1878  1.1892  0.3918
   lcl|BK006935.2_cds_DAA06921.1_4 1.1082  0.9759  0.3646
   lcl|BK006935.2_cds_DAA06922.1_5 1.0679  1.0     0.3447
   lcl|BK006935.2_cds_DAA06923.1_6 1.1587  1.1892  0.3668
   lcl|BK006935.2_cds_DAA06924.2_7 1.1429  1.1892  0.5033
   lcl|BK006935.2_cds_DAA06925.1_8 1.2349  1.2315  0.4547
   lcl|BK006935.2_cds_DAA06926.1_9 1.1726  1.2315  0.3764
   ...                             ...     ...     ...

Here, only the first nine lines of the output are shown. The first column is the
gene identifier, the second, third, and fourth columns are the mean, median, and
standard deviation of RSCU values observed in the gene. Our analyses demonstrate
that the mean gw-RSCU performs best to examine codon optimization. These metrics
can also be used to examine gene-wise codon usage biases.

|

.. |br| raw:: html

  <br/>
