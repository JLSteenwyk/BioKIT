.. _change_log:


Change log
==========

^^^^^

Major changes to BioKIT are summarized here.

*0.1.0*: Functions that remove adapters from FASTQ files (trim_se_adapters_fastq
and trim_pe_adapters_fastq) are now part of BioKIT. Note, these identify
exact matches of adapters. 

Parsimony informative sites, constant sites, and variable sites functions
now have a verbose option that allows users to examine the characterization
of each site in an alignment.

The name gw-RSCU has been shortened to gRSCU.

*0.0.9*: Functions that look at codons (e.g., RSCU and gw-RSCU) now can account for ambiguous codons.
For example, codons that have ambiguous characters like the codon "CNN." These codons
are skipped during analysis of RSCU and gw-RSCU.

*1.0.1*: Added "X" as a gap character during alignment recoding

*1.1.0*: Added Dayhoff-9, -12, -15, and -18 recoding schemes 

*1.1.3*: Fixed bug in RSCU calculations
