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

Alignment length 
################
Function names: alignment_length; aln_len |br|
Command line interface: bk_alignment_length; bk_aln_len

Calculate the length of an alignment. 

.. code-block:: shell

	biokit alignment_length <fasta>

Options: |br|
*<alignment>*: first argument after function name should be an alignment file 

|

Alignment recoding 
##################
Function names: alignment_recoding; aln_recoding; recode |br|
Command line interface: bk_alignment_recoding; bk_aln_recoding; bk_recode

Recode alignments using reduced character states.

Alignments can be recoded using established or
custom recoding schemes. Recoding schemes are
specified using the -c/--code argument. Custom
recoding schemes can be used and should be formatted
as a two column file wherein the first column is the
recoded character and the second column is the character
in the alignment.

.. code-block:: shell

	biokit alignment_recoding <fasta> [-c/--code <code>]

Codes for which recoding scheme to use: |br|
**RY-nucleotide** |br|
R = purines (i.e., A and G) |br|
Y = pyrimidines (i.e., T and C)

**Dayhoff-6** |br|
0 = A, G, P, S, and T |br|
1 = D, E, N, and Q |br|
2 = H, K, and R |br|
3 = I, L, M, and V |br|
4 = F, W, and Y |br|
5 = C

**SandR-6** |br|
0 = A, P, S, and T |br|
1 = D, E, N, and G |br|
2 = Q, K, and R |br|
3 = M, I, V, and L |br|
4 = W and C |br|
5 = F, Y, and H |br|

**KGB-6** |br|
0 = A, G, P, and S |br|
1 = D, E, N, Q, H, K, R, and T |br|
2 = M, I, and L |br|
3 = W |br|
4 = F and Y |br|
5 = C and V |br|

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-c/\\-\\-code*: argument to specify the recoding scheme to use 

|

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

Constant sites
##############
Function names: constant_sites; con_sites |br|
Command line interface: bk_constant_sites; bk_con_sites

Calculate the number of constant sites in an alignment.

Constant sites are defined as a site in an
alignment with the same nucleotide or amino
acid sequence (excluding gaps) among all taxa.

.. code-block:: shell

	biokit constant_sites <fasta>

Options: |br|
*<alignment>*: first argument after function name should be an alignment file 

|

Parsimony informative sites
###########################
Function names: parsimony_informative_sites; pi_sites; pis |br|
Command line interface: bk_parsimony_informative_sites; bk_pi_sites; bk_pis

Calculate the number of parsimony informative
sites in an alignment.

Parsimony informative sites are defined as a
site in an alignment with at least two nucleotides
or amino acids that occur at least twice.

.. code-block:: shell

	biokit parsimony_informative_sites <fasta>

Options: |br|
*<alignment>*: first argument after function name should be an alignment file 

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

Variable sites
##############
Function names: variable_sites; var_sites; vs |br|
Command line interface: bk_variable_sites; bk_var_sites; bk_vs

Calculate the number of variable sites in an
alignment.

Variable sites are defined as a site in an
alignment with at least two nucleotide or amino
acid characters among all taxa.

.. code-block:: shell

	biokit variable_sites <fasta>

Options: |br|
*<alignment>*: first argument after function name should be an alignment file 

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

Gene-wise relative synonymous codon usage
#########################################
Function names: gene_wise_relative_synonymous_codon_usage; gene_wise_rscu; gw_rscu |br|
Command line interface: bk_gene_wise_relative_synonymous_codon_usage; bk_gene_wise_rscu; bk_gw_rscu

Calculate gene-wise relative synonymous codon usage (gw-RSCU).

Codon usage bias examines biases for codon usage of
a particular gene. We adapted RSCU to be applied to
individual genes rather than only codons. Specifically,
gw-RSCU is the mean (or median) RSCU value observed
in a particular gene. This provides insight into how
codon usage bias influences codon usage for a particular
gene. This function also outputs the standard deviation
of RSCU values for a given gene.

The output is col 1: the gene identifier; col 2: the
gw-RSCU based on the mean RSCU value observed in a gene;
col 3: the gw-RSCU based on the median RSCU value observed
in a gene; and the col 4: the standard deviation of
RSCU values observed in a gene.

Custom genetic codes can be used as input and should
be formatted with the codon in first column and the 
resulting amino acid in the second column.

.. code-block:: shell

   biokit gene_wise_relative_synonymous_codon_usage <fasta> [-tt/--translation_table <code>]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-tt/\\-\\-translation_table*: optional argument of the code for the translation
table to be used. Default: 1, which is the standard code.

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

   biokit relative_synonymous_codon_usage <fasta> [-tt/--translation_table <code>]

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

N50 is the sequence length of the shortest contig at 50% of the genome size.

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

   biokit n90 <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file

|

Number of large scaffolds
#########################
Function names: number_of_large_scaffolds; num_of_lrg_scaffolds; number_of_large_contigs; num_of_lrg_cont |br|
Command line interface: bk_number_of_large_scaffolds; bk_num_of_lrg_scaffolds; bk_number_of_large_contigs; bk_num_of_lrg_cont

Calculate number and total sequence length of
large scaffolds. Each value is represented as
column 1 and column 2 in the output, respectively.

.. code-block:: shell

   biokit number_of_large_scaffolds <fasta> [-t/--threshold <int>]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file
*-t/\\-\\-threshold*: threshold for what is considered
a large scaffold. Only scaffolds with a length greater than this
value will be counted. Default: 500

|

Number of scaffolds
###################
Function names: number_of_scaffolds; num_of_scaffolds; number_of_contigs; num_of_cont |br|
Command line interface: bk_number_of_scaffolds; bk_num_of_scaffolds; bk_number_of_contigs; bk_num_of_cont

Calculate the number of scaffolds or entries
in a FASTA file. In this way, a user can also 
determine the number of predicted genes in a 
coding sequence or protein FASTA file with this
function.

.. code-block:: shell

   biokit number_of_scaffolds <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file

|

Sum of scaffold lengths
#######################
Function names: sum_of_scaffold_lengths; sum_of_contig_lengths |br|
Command line interface: bk_sum_of_scaffold_lengths; bk_sum_of_contig_lengths

Determine the sum of scaffold lengths. 
                
The intended use of this function is to determine
the length of a genome assembly, but can also be
used, for example, to determine the sum length
of all coding sequences.

.. code-block:: shell

   biokit sum_of_scaffold_lengths <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file

|

Sequence summary and processing functions
-----------------------------------------

^^^^^

Character frequency
###################
Function names: character_frequency; char_freq |br|
Command line interface: bk_character_frequency; bk_char_freq

Calculate the frequency of characters in a FASTA file.

This can be used to determine the frequency of A, T, C, and G
in a genome or the frequency of amino acids in a proteome.

.. code-block:: shell

   biokit character_frequency <fasta> [-v/--verbose]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|

|

Get FASTA entry (faidx)
#######################
Function names: faidx; get_entry; ge |br|
Command line interface: bk_faidx; bk_get_entry; bk_ge

Extracts sequence entry from fasta file.

This function works similarly to the faidx function 
in samtools, but does not requiring an indexing the
sequence file.

.. code-block:: shell

   biokit faidx <fasta> [-v/--verbose]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-e/\\-\\-entry*: entry name to be extracted from the inputted fasta file

|

File format converter
#####################
Function names: file_format_converter; format_converter; ffc |br|
Command line interface: bk_file_format_converter; bk_format_converter; bk_ffc

Converts a multiple sequence file from one format to another.

Acceptable file formats include FASTA, Clustal, MAF, Mauve,
Phylip, Phylip-sequential, Phylip-relaxed, and Stockholm.
Input and output file formats are specified with the
\\-\\-input_file_format and \\-\\-output_file_format arguments;
input and output files are specified with the \\-\\-input_file
and \\-\\-output_file arguments.

.. code-block:: shell

   biokit file_format_converter -i/--input_file <input_file> -iff/--input_file_format <input_file_format>  -o/--output_file <output_file> -off/--output_file_format <output_file_format>

Options: |br|
*-i/\\-\\-input_file*: input file name 
*-iff/\\-\\-input_file_format*: input file format
*-o/\\-\\-output_file*: output file name
*-off/\\-\\-output_file_format*: output file format

|

Multiple line to single line FASTA
##################################
Function names: multiple_line_to_single_line_fasta; ml2sl |br|
Command line interface: bk_multiple_line_to_single_line_fasta; bk_ml2sl

Converts FASTA files with multiple lines
per sequence to a FASTA file with the sequence
represented on one line.

.. code-block:: shell

   biokit multiple_line_to_single_line_fasta <fasta> [-o/--output <output_file>]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-o/\\-\\-output*: optional argument to name the output file

|

Remove FASTA entry
##################
Function names: remove_fasta_entry |br|
Command line interface: bk_remove_fasta_entry

Remove FASTA entry from multi-FASTA file.

Output will have the suffix "pruned.fa" unless
the user specifies a different output file name. 

.. code-block:: shell

   biokit remove_fasta_entry <fasta> -e/--entry <entry> [-o/--output <output_file>]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-e/\\-\\-entry*: entry name to be removed from the inputted fasta file |br|
*-o/\\-\\-output*: optional argument to write the renamed fasta file to.
Default output has the same name as the input file with the suffix "pruned.fa"
added to it.

|

Remove short sequences
######################
Function names: remove_short_sequences; remove_short_seqs |br|
Command line interface: bk_remove_short_sequences; bk_remove_short_seqs

Remove short sequences from a multi-FASTA file.

Short sequences are defined as having a length
less than 500. Users can specify their own threshold.
All sequences greater than the threshold will be
kept in the resulting file.

Output will have the suffix "long_seqs.fa" unless
the user specifies a different output file name.

.. code-block:: shell

   biokit remove_short_sequences <fasta> -t/--threshold <threshold> [-o/--output <output_file>]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-t/\\-\\-threshold*: threshold for short sequences. Sequences greater 
than this value will be kept |br|
*-o/\\-\\-output*: optional argument to write the renamed fasta file to.
Default output has the same name as the input file with the suffix "long_seqs.fa"
added to it.

|

Rename FASTA entries
####################
Function names: rename_fasta_entries; rename_fasta |br|
Command line interface: bk_rename_fasta_entries; bk_rename_fasta

Renames fasta entries.

Renaming fasta entries will follow the scheme of a tab-delimited
file wherein the first column is the current fasta entry name and
the second column is the new fasta entry name in the resulting 
output alignment. 

.. code-block:: shell

   rename_fasta_entries <fasta> -i/--idmap <idmap> [-o/--output <output_file>]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-i/\\-\\-idmap*: identifier map of current FASTA names (col1) and desired FASTA names (col2) |br|
*-o/\\-\\-output*: optional argument to name the output file

|

Reorder by sequence length
##########################
Function names: reorder_by_sequence_length; reorder_by_seq_len |br|
Command line interface: bk_reorder_by_sequence_length; bk_reorder_by_seq_len

Reorder FASTA file entries from the longest entry
to the shortest entry. 

.. code-block:: shell

   biokit reorder_by_sequence_length <fasta> [-o/--output <output_file>]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-o/\\-\\-output*: optional argument to write the reordered fasta file to.
Default output has the same name as the input file with the suffix
".reordered.fa" added to it.

|

Sequence complement
###################
Function names: sequence_complement; seq_comp |br|
Command line interface: bk_sequence_complement; bk_seq_comp

Generates the sequence complement for all entries
in a multi-FASTA file. To generate a reverse sequence
complement, add the -r/--reverse argument.

.. code-block:: shell

   biokit sequence_complement <fasta> [-r/--reverse]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-r/\\-\\-reverse*: if used, the reverse complement sequence will be generated

|

Sequence length
###############
Function names: sequence_length; seq_len |br|
Command line interface: bk_sequence_length; bk_seq_len

Calculate sequence length of each FASTA entry.

.. code-block:: shell

   biokit sequence_length <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|

|

Single line to multiple line fasta
##################################
Function names: single_line_to_multiple_line_fasta; sl2ml |br|
Command line interface: bk_single_line_to_multiple_line_fasta; bk_sl2ml

Calculate sequence length of each FASTA entry.

.. code-block:: shell

   biokit single_line_to_multiple_line_fasta <fasta>

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|

|

.. |br| raw:: html

  <br/>
