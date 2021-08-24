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

   phykit aln_summary -h

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

.. |br| raw:: html

  <br/>
