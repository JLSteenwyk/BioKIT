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

*1.1.4*: Added structured output support (tsv/json/yaml) to additional commands and
expanded integration coverage with output format consistency and CLI output
contract tests.

*1.1.5*: Updated dependency pinning to fix installation failures and moved supported
Python versions to 3.11, 3.12, and 3.13 only.

*1.1.19*: Added ``neutrality_plot`` command for GC12 vs. GC3 regression
analysis of coding sequences. The regression slope is a classical codon
usage diagnostic — a slope near 1 indicates that codon usage is largely
mutation-driven, while a slope near 0 indicates that GC12 is constrained
by selection while GC3 drifts. Reports slope, intercept, r-squared, and
n; with ``-v/--verbose`` includes per-gene GC12 and GC3 values; with
``-p/--plot`` saves a scatter plot with the regression line as PNG.

*1.1.18*: Added ``sample_sequences`` (alias: ``sample_seqs``) command to
randomly draw N sequences (``-n/--number``) or a percentage
(``-p/--percent``, default 10%) from a FASTA file without replacement.
Supports reproducible sampling via ``-s/--seed`` and optional output
file via ``-o/--output``.

*1.1.17*: Added ``find_orfs`` (alias: ``orfs``) command to identify all open
reading frames above a minimum length in nucleotide sequences, searching
all 6 reading frames (3 forward, 3 reverse). Reports id, frame, start,
stop, and length in nucleotides and amino acids. With ``--extract``,
outputs the ORF sequences as FASTA; with ``--protein``, outputs protein
translations instead of nucleotides. Supports custom genetic codes via
``-tt/--translation_table``.

*1.1.16*: Added ``kmer_frequency`` (alias: ``kmer_freq``) command to count
all k-mers of a given size in sequences. Useful for genomic signatures,
metagenomics binning prep, and alignment-free sequence comparison.
Supports ``--canonical`` to collapse each k-mer with its reverse complement
and ``-v/--verbose`` for per-sequence counts. K-mers containing non-ACGT
characters are skipped.

*1.1.15*: Added ``homopolymer_runs`` (alias: ``homopolymer``) command to
find the longest homopolymer run per sequence in a FASTA file. Reports
length, base, and 1-based start position. Relevant for nanopore/PacBio
QC where homopolymer errors are the dominant error mode. Use
``--per-base`` to report the longest run for each of A, C, G, T.

*1.1.14*: Added ``fastq_to_fasta`` (alias: ``fq2fa``) command to convert
FASTQ files to FASTA by stripping quality scores. Supports stdin input
via ``-`` and optional output file via ``-o/--output``.

*1.1.13*: Added ``genbank_to_fasta`` (alias: ``gb2fa``) command to extract
sequences from GenBank flat files. Optionally filters by feature type
(e.g., CDS, rRNA, tRNA, gene) via ``-t/--feature_type`` and can output
protein translations for CDS features via ``--translate``.

*1.1.12*: Added ``assembly_curve`` (alias: ``asm_curve``) command to output
cumulative length vs. contig rank (sorted by descending length). Supports
optional plotting via ``-p/--plot`` to generate a PNG assembly curve figure.

*1.1.11*: Added ``melting_temperature`` (alias: ``tm``) command to calculate melting
temperature of nucleotide sequences using nearest-neighbor thermodynamics. Supports
custom Na+ and oligo concentrations.

*1.1.10*: Added ``shuffle_sequences`` (alias: ``shuffle_seqs``) command to randomly
shuffle nucleotide or amino acid order within each sequence while preserving
composition. Supports reproducible shuffling via ``-s/--seed``.

*1.1.9*: Added ``restriction_sites`` (alias: ``re_sites``) command to find restriction
enzyme recognition sites in sequences. Reports cut positions and fragment sizes
for one or more enzymes per sequence. Uses BioPython's Restriction module.

*1.1.8*: Added ``gc_content_four_fold_degenerate_sites`` (alias: ``gc4``) command to
calculate GC content at four-fold degenerate sites in coding sequences. Supports
custom translation tables via ``-tt`` and per-sequence output via ``-v``.

*1.1.7*: Added ``protein_charge`` (alias: ``prot_charge``) command to calculate the
net charge of protein sequences at a given pH (default 7.0). Uses BioPython's
ProteinAnalysis for the calculation. Supports tsv, json, and yaml output.

*1.1.6*: Added ``fasta_deduplication`` (alias: ``dedup``) command to remove duplicate
sequences from multi-FASTA files. Sequences are compared case-insensitively
and the first occurrence of each unique sequence is kept.
