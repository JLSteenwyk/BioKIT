.. _usage:

Usage
=====

BioKIT helps process and analyze data for bioinformatics research.

Generally, all functions are designed to help understand the contents of
various files including those representing alignments, coding sequences, fastq
files, and genomes.

|

General usage
-------------

^^^^^

Calling functions
#################

.. code-block:: shell

   biokit <command> [optional command arguments]

Command specific help messages can be viewed by adding a 
-h/\\-\\-help argument after the command. For example, to see the
to see the help message for the command 'alignment_summary', execute:

.. code-block:: shell

   biokit alignment_summary -h 
   # or
   biokit alignment_summary --help

|

Function aliases
################

Each function comes with aliases to save the user some
key strokes. For example, to get the help message for the 'alignment_summary'
function, you can type:

.. code-block:: shell

   biokit aln_summary -h

|

Command line interfaces
#######################

All functions (including aliases) can also be executed using
a command-line interface that starts with *bk_*. For example, instead of typing
the previous command to get the help message of the alignment_summary function,
you can type:

.. code-block:: shell

   bk_alignment_summary -h
   # or
   bk_aln_summary -h


All possible function names are specified at the top of each function section. 

|

Alignment-based functions
-------------------------

^^^^^

Alignment summary 
#################
Function names: alignment_summary; aln_summary |br|
Command line interface: bk_alignment_summary; bk_aln_summary

Summary statistics for an alignment. Reported
statistics include alignment length, number of taxa,
number of parsimony sites, number of variable sites,
number of constant sites, frequency of each character
(including gaps, which are considered to be '-' or '?'). 

.. code-block:: shell

	biokit alignment_summary <fasta>

Options: |br|
*<alignment>*: first argument after function name should be an alignment file 

|

Consensus sequence
##################
Function names: consensus_sequence; con_seq |br|
Command line interface: bk_consensus_sequence; bk_con_seq

Generates a consequence from a multiple sequence alignment file in FASTA format.

.. code-block:: shell

	biokit consensus_sequence <fasta> -t/--threshold <threshold> -ac/--ambiguous_character <ambiguous character>

Options: |br|
*<fasta>*: first argument after function name should be an alignment fasta file 
*<threshold>*: threshold for how common a residue must be to be represented  
*<ambiguous character>*: the ambiguity character to use. Default is 'N'

|

Position specific score matrix
##############################
Function names: position_specific_score_matrix; pssm |br|
Command line interface: bk_position_specific_score_matrix; bk_pssm

Generates a position specific score matrix for an alignment.

.. code-block:: shell

	biokit position_specific_score_matrix <fasta> [-ac/--ambiguous_character <ambiguous character>]

Options: |br|
*<fasta>*: first argument after function name should be an alignment fasta file 
*<ambiguous character>*: the ambiguity character to use. Default is 'N'

|

Coding sequence-based functions
-------------------------------

^^^^^

GC content first codon position
###############################
Function names: gc_content_first_position; gc1 |br|
Command line interface: bk_gc_content_first_position; bk_gc1

Calculate GC content of the first codon position.
The input must be the coding sequence of a gene or
genes. All genes are assumed to have sequence lengths
divisible by three.

.. code-block:: shell

   biokit gc_content_first_position <fasta> [-v/--verbose]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-v/\\-\\-verbose*: optional argument to print the GC content of each fasta entry

|

GC content second codon position
###############################
Function names: gc_content_second_position; gc2 |br|
Command line interface: bk_gc_content_second_position; bk_gc2

Calculate GC content of the second codon position.
The input must be the coding sequence of a gene or
genes. All genes are assumed to have sequence lengths
divisible by three.

.. code-block:: shell

   biokit gc_content_second_position <fasta> [-v/--verbose]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-v/\\-\\-verbose*: optional argument to print the GC content of each fasta entry

|

GC content third codon position
###############################
Function names: gc_content_third_position; gc3 |br|
Command line interface: bk_gc_content_third_position; bk_gc3

Calculate GC content of the third codon position.
The input must be the coding sequence of a gene or
genes. All genes are assumed to have sequence lengths
divisible by three.

.. code-block:: shell

   biokit gc_content_third_position <fasta> [-v/--verbose]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-v/\\-\\-verbose*: optional argument to print the GC content of each fasta entry

|

Relative synonymous codon usage
###############################
Function names: relative_synonymous_codon_usage; rscu |br|
Command line interface: bk_relative_synonymous_codon_usage; bk_rscu

Calculate relative synonymous codon usage.

Relative synonymous codon usage is the ratio
of the observed frequency of codons over the
expected frequency given that all the synonymous
codons for the same amino acids are used equally.

Custom genetic codes can be used as input and should
be formatted with the codon in first column and the 
resulting amino acid in the second column.

.. code-block:: shell

   biokit relative_synonymous_codon_usage <fasta> [-v/--verbose]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-tt/\\-\\-translation_table*: optional argument of the code for the translation
table to be used. Default: 1, which is the standard code.

|

Translate sequence
##################
Function names: translate_sequence; translate_seq; trans_seq |br|
Command line interface: bk_translate_sequence; bk_translate_seq; bk_trans_seq

Translates coding sequences to amino acid
sequences. Sequences can be translated using
diverse genetic codes. For codons that can
encode two amino acids (e.g., TAG encodes
Glu or STOP in the Blastocrithidia Nuclear Code),
the standard genetic code is used.

Custom genetic codes can be used as input and should
be formatted with the codon in first column and the 
resulting amino acid in the second column.

.. code-block:: shell

   biokit translate_sequence <fasta> [-tt/--translation_table <code> -o/--output <output_file>]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-tt/\\-\\-translation_table*: optional argument of the code for the translation
table to be used. Default: 1, which is the standard code.
*-o/\\-\\-output*: optional argument to write the translated fasta file to.
Default output has the same name as the input file with the suffix ".translated.fa" added
to it.

|

FASTQ file functions
--------------------

^^^^^

FASTQ read lengths
##################
Function names: fastq_read_lengths; fastq_read_lens |br|
Command line interface: bk_fastq_read_lengths; bk_fastq_read_lens

Determine lengths of FASTQ reads.
                
Using default arguments, the average and
standard deviation of read lengths in a
FASTQ file will be reported. To obtain
the lengths of all FASTQ reads, use the
verbose option.

.. code-block:: shell

   biokit fastq_read_lengths <fasta> [-tt/--translation_table <code> -o/--output <output_file>]

Options: |br|
*<fastq>*: first argument after function name should be a FASTQ file |br|
*-v/\\-\\-verbose*: print length of each FASTQ read

|

Subset PE FASTQ reads
#####################
Function names: subset_pe_fastq_reads; subset_pe_fastq |br|
Command line interface: bk_subset_pe_fastq_reads; bk_subset_pe_fastq

Subset paired-end FASTQ data.

Subsetting FASTQ data may be helpful for 
running test scripts or achieving equal 
coverage between samples. A percentage of
total reads in paired-end FASTQ data can
be obtained with this function. Random
subsamples are obtained using seeds for
reproducibility. If no seed is specified,
a seed is generated based off of the date
and time. During subsampling, paired-end
information is maintained.

Files are outputed with the suffix "_subset.fq"

.. code-block:: shell

   biokit subset_pe_fastq_reads <fastq1> <fastq2> [-p/--percent <percent> -s/--seed <seed> -o/--output_file <output_file>]

Options: |br|
*<fastq1>*: first argument after function name should be the name of one of the fastq files |br|
*<fastq2>*: first argument after function name should be the name of the other fastq file |br|
*-p/\\-\\-percent*: percentage of reads to maintain in subsetted data |br|
*-s/\\-\\-seed*: seed for random sampling

|

Subset SE FASTQ reads
#####################
Function names: subset_se_fastq_reads; subset_se_fastq |br|
Command line interface: bk_subset_se_fastq_reads; bk_subset_se_fastq

Subset single-end FASTQ data.

Subsetting FASTQ data may be helpful for 
running test scripts or achieving equal 
coverage between samples. A percentage of
total reads in single-end FASTQ data can
be obtained with this function. Random
subsamples are obtained using seeds for
reproducibility. If no seed is specified,
a seed is generated based off of the date
and time.

Output files will have the suffix "_subset.fq"

.. code-block:: shell

   biokit subset_se_fastq_reads <fastq> [-p/--percent <percent> -s/--seed <seed> -o/--output_file <output_file>]

Options: |br|
*<fastq>*: first argument after function name should be the name of one of the fastq files |br|
*-p/\\-\\-percent*: percentage of reads to maintain in subsetted data |br|
*-s/\\-\\-seed*: seed for random sampling |br|
*-o/\\-\\-output_file*: specify the name of the output file

|

Trim PE FASTQ reads
###################
Function names: trim_pe_fastq_reads; trim_pe_fastq |br|
Command line interface: bk_trim_pe_fastq_reads; bk_trim_pe_fastq

Quality trim paired-end FastQ data.

FASTQ data will be trimmed according to
quality score and length of the reads.
Users can specify quality and length
thresholds. Paired reads that are 
maintained and saved to files with the
suffix "_paired_trimmed.fq." Single
reads that passed quality thresholds are
saved to files with the suffix 
"_unpaired_trimmed.fq."

.. code-block:: shell

   biokit trim_pe_fastq_reads <fastq1> <fastq2> [-m/--minimum 20 -l/--length 20]

Options: |br|
*<fastq1>*: first argument after function name should be the name of one of the fastq files |br|
*<fastq2>*: first argument after function name should be the name of the other fastq file |br|
*-m/\\-\\-minimum*: minimum quality of read to be kept. Default: 20 |br|
*-l/\\-\\-length*: minimum length of read to be kept. Default: 20

|

Trim SE FASTQ reads
###################
Function names: trim_se_fastq_reads; trim_se_fastq |br|
Command line interface: bk_trim_se_fastq_reads; bk_trim_se_fastq

Quality trim single-end FastQ data.

FASTQ data will be trimmed according to
quality score and length of the reads.
Users can specify quality and length
thresholds. Output file has the suffix 
"_trimmed.fq" or can be named by the user 
with the output_file argument.

.. code-block:: shell

   biokit trim_se_fastq_reads <fastq> [-m/--minimum 20 -l/--length 20]

Options: |br|
*<fastq>*: first argument after function name should be the fastq file |br|
*-m/\\-\\-minimum*: minimum quality of read to be kept. Default: 20 |br|
*-l/\\-\\-length*: minimum length of read to be kept. Default: 20 |br|
*-o/\\-\\-output_file*: name of the output file of trimmed reads

|

Genome functions
----------------

^^^^^

GC content
##########
Function names: gc_content; gc |br|
Command line interface: bk_gc_content; bk_gc

Calculate GC content of a fasta file.

GC content is the fraction of bases that are
either guanines or cytosines. To obtain
GC content per FASTA entry, use the verbose
option.

.. code-block:: shell

   biokit gc_content <fasta> [-v/--verbose]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-v/\\-\\-verbose*: optional argument to print the GC content of each fasta entry

|

Genome assembly metrics
#######################
Function names: genome_assembly_metrics; assembly_metrics |br|
Command line interface: bk_genome_assembly_metrics; bk_assembly_metrics

Calculate L50, L90, N50, N90, GC content, assembly size,
number of scaffolds, number and sum length
of large scaffolds, frequency of A, T, C, and G.

L50: The smallest number of contigs whose length sum makes up half of the genome size. |br|
L90: The smallest number of contigs whose length sum makes up 90% of the genome size. |br|
N50: The sequence length of the shortest contig at half of the genome size. |br|
N90: The sequence length of the shortest contig at 90% of the genome size. |br|
GC content: The fraction of bases that are either guanines or cytosines. |br|
Assembly size: The sum length of all contigs in an assembly. |br|
Number of scaffolds: The total number of scaffolds in an assembly. |br|
Number of large scaffolds: The total number of scaffolds that are greater than the threshold for small scaffolds. |br|
Sum length of large scaffolds: The sum length of all large scaffolds. |br|
Frequency of A: The number of occurences of A corrected by assembly size. |br|
Frequency of T: The number of occurences of T corrected by assembly size. |br|
Frequency of C: The number of occurences of C corrected by assembly size. |br|
Frequency of G: The number of occurences of G corrected by assembly size.

.. code-block:: shell

   biokit genome_assembly_metrics <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-t/\\-\\-threshold*: threshold for what is considered a large scaffold.
Only scaffolds with a length greater than this value will be counted. Default: 500

|

L50
###
Function names: l50 |br|
Command line interface: bk_l50

Calculates L50 for a genome assembly.

L50 is the smallest number of contigs whose length sum
makes up half of the genome size.

.. code-block:: shell

   biokit l50 <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file

|

L90
###
Function names: l90 |br|
Command line interface: bk_l90

Calculates L90 for a genome assembly.

L90 is the smallest number of contigs whose length sum
makes up 90% of the genome size.

.. code-block:: shell

   biokit l90 <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file

|

Longest scaffold
################
Function names: longest_scaffold; longest_scaff; longest_contig; longest_cont |br|
Command line interface: bk_longest_scaffold; bk_longest_scaff; bk_longest_contig; bk_longest_cont

Determine the length of the longest scaffold in a genome assembly.

.. code-block:: shell

   biokit longest_scaffold <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file

|

N50
###
Function names: n50 |br|
Command line interface: bk_n50

Calculates N50 for a genome assembly.

N50 is the smallest number of contigs whose length sum
makes up half of the genome size.

.. code-block:: shell

   biokit n50 <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file

|

N90
###
Function names: n90 |br|
Command line interface: bk_n90

Calculates N90 for a genome assembly.

N90 is the sequence length of the shortest contig at 90% of the genome size.

.. code-block:: shell

   biokit l90 <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file

|

.. |br| raw:: html

  <br/>
