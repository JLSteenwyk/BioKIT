#!/usr/bin/env python

import logging
import sys
import textwrap
from .version import __version__

from argparse import (
    ArgumentParser,
    RawTextHelpFormatter,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .services.text import (
    CharacterFrequency,
    ConsensusSequence,
    Faidx,
    FileFormatConverter,
    GCContent,
    GCContentThirdPosition,
    GenomeAssemblyMetrics,
    L50,
    L90,
    LongestScaffold,
    N50,
    N90,
    NumberOfScaffolds,
    NumberOfLargeScaffolds,
    PositionSpecificScoreMatrix,
    RenameFastaEntries,
    ReorderBySequenceLength,
    SequenceComplement,
    SequenceLength,
    SumOfScaffoldLengths,
    TranslateSequence,
    TrimSEFastQ,
)

from .services.tree import (
    TipLabels,
)

from .helpers.boolean_argument_parsing import str2bool


logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

help_header = f"""
                 ____  _       _  _______ _______ 
                |  _ \(_)     | |/ /_   _|__   __|
                | |_) |_  ___ | ' /  | |    | |   
                |  _ <| |/ _ \|  <   | |    | |   
                | |_) | | (_) | . \ _| |_   | |   
                |____/|_|\___/|_|\_\_____|  |_|   
                            
                Version: {__version__}
                Citation: Steenwyk et al. 2021, CITATION INFORMATION
                
"""

class Biokit(object):
    help_header = f"""
                 ____  _       _  _______ _______ 
                |  _ \(_)     | |/ /_   _|__   __|
                | |_) |_  ___ | ' /  | |    | |   
                |  _ <| |/ _ \|  <   | |    | |   
                | |_) | | (_) | . \ _| |_   | |   
                |____/|_|\___/|_|\_\_____|  |_|   
                            
                Version: {__version__}
                Citation: Steenwyk et al. 2021, CITATION INFORMATION

    """
    
    def __init__(self):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                BioKIT is a broadly applicable command-line toolkit for bioinformatics research.

                Usage: biokit <command> [optional command arguments]

                Command specific help messages can be viewed by adding a 
                -h/--help argument after the command. For example, to see the
                to see the help message for the command 'get_entry', execute
                "biokit get_entry -h" or "biokit get_entry --help".

                Lastly, each function comes with aliases to save the user some
                key strokes. For example, to get the help message for the 'get_entry'
                function, you can type "biokit ge -h". All aliases are specified
                in parentheses after the long form of the function name. 

                Text sequence file-based commands
                =================================
                consensus_sequence (alias: con_len)
                    - create a consensus sequence from an alignment
                faidx (alias: get_entry; ge)
                    - extract query fasta entry from multi-fasta file
                file_format_converter (alias: format_converter; ffc)
                    - convert a multiple sequence file from one format
                      to another
                gc_content (alias: gc)
                    - calculate the GC content of a FASTA file
                gc_content_third_position (alias: gc3)
                    - calculate the GC content of the third position
                      among coding sequences
                l50
                    - calculate the L50 of a genome assembly
                l90
                    - calcualte the L90 of a genome assembly
                longest_scaffold (alias: longest_scaff, longest_contig, longest_cont)
                    - determine the length of the longest
                      scaffold of a genome assembly
                n50
                    - calculate the N50 of a genome assembly
                n90
                    - calculate the N90 of a genome assembly
                number_of_scaffolds (alias: num_of_scaffolds, number_of_contigs, num_of_cont)
                    - calculate the number of scaffolds in a
                      genome assembly
                number_of_large_scaffolds (alias: num_of_lrg_scaffolds, number_of_large_contigs, num_of_lrg_cont)
                    - calculate the number of large scaffolds
                position_specific_score_matrix (alias: pssm)
                    - create a position specific score matrix for an alignment
                rename_fasta_entries (alias: rename_fasta)
                    - rename entries in a FASTA file
                reorder_by_sequence_length (alias: reorder_by_seq_len)
                    - reorder sequences from longest to shortest in a FASTA file
                sequence_complement (alias: seq_comp)
                    - generate the complementary sequence for an alignment 
                sequence_length (alias: seq_len)
                    - calculate the length of each FASTA entry
                translate_sequence (alias: translate_seq, trans_seq) 
                    - translate coding sequences to amino acids

                Tree-based commands
                ===================
                tip_labels (alias: tree_labels; labels; tl)
                    - print leaf names in a phylogeny
                
                """
            ),
        )
        parser.add_argument("command", help=SUPPRESS)
        args = parser.parse_args(sys.argv[1:2])

        # if command is part of the possible commands (i.e., the long form
        # commands, run). Otherwise, assume it is an alias and look to the
        # run_alias function
        try:
            if hasattr(self, args.command):
                getattr(self, args.command)(sys.argv[2:])
            else:
                self.run_alias(args.command, sys.argv[2:])
        except NameError:
            sys.exit()

    ## Aliases
    def run_alias(self, command, argv):
        # version
        if command in ['v']:
            return self.version()
        # Text aliases
        elif command in ['con_len']:
            return self.consensus_sequence(argv)
        elif command in ['get_entry', 'ge']:
            return self.faidx(argv)
        elif command in ['format_converter', 'ffc']:
            return self.file_format_converter(argv)
        elif command in ['longest_scaff', 'longest_contig', 'longest_cont']:
            return self.longest_scaffold(argv)
        elif command in ['num_of_scaffolds', 'number_of_contigs', 'num_of_cont']:
            return self.number_of_scaffolds(argv)
        elif command in [
            'number_of_large_scaffolds',
            'num_of_lrg_scaffolds', 
            'number_of_large_contigs',
            'num_of_lrg_cont'
        ]:
            return self.number_of_large_scaffolds(argv)
        elif command in ['pssm']:
            return self.position_specific_score_matrix(argv)
        elif command in ['reorder_by_seq_len']:
            return self.reorder_by_sequence_length(argv)
        elif command in ['seq_comp']:
            return self.sequence_complement(argv)
        elif command in ['seq_len']:
            return self.sequence_length(argv)
        # Tree aliases
        elif command in ['labels', 'tree_labels', 'tl']:
            return self.tip_labels(argv)
        else:
            print("Invalid command option. See help for a complete list of commands and aliases.")
            parser.print_help()
            sys.exit(1)

    ## print version
    def version(self):
        print(textwrap.dedent(
            f"""\
            {self.help_header}
            """
        ))

    ## Alignment functions
    @staticmethod
    def character_frequency(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate the frequency of characters in a FASTA file.
                
                Aliases:
                  character_frequency, char_freq
                Command line interfaces: 
                  bk_character_frequency, bk_char_freq

                Usage:
                biokit character_frequency <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        CharacterFrequency(args).run()

    @staticmethod
    def consensus_sequence(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Generates a consequence from a multiple sequence alignment
                file in FASTA format.
                
                Aliases:
                  consensus_sequence, con_len
                Command line interfaces: 
                  bk_consensus_sequence, bk_con_len

                Usage:
                biokit consensus_sequence <fasta> -t/--threshold <threshold>
                -ac/--ambiguous_character <ambiguous character>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -t/--threshold              threshold for how common
                                            a residue must be to be
                                            represented
                
                -ac/--ambiguous_character   the ambiguity character to
                                            use. Default is 'N'
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-t","--threshold", type=str, help=SUPPRESS)
        parser.add_argument("-ac","--ambiguous_character", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        ConsensusSequence(args).run()

    @staticmethod
    def faidx(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Extracts sequence entry from fasta file.

                This function works similarly to the faidx function 
                in samtools, but does not requiring an indexing function.

                Aliases:
                  faidx, get_entry, ge
                Command line interfaces: 
                  bk_faidx, bk_get_entry, bk_ge

                Usage:
                biokit faidx <fasta> -e/--entry <fasta entry>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be a
                                            query fasta file

                -e/--entry                  entry name to be extracted
                                            from the inputted fasta file
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-e","--entry", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        Faidx(args).run()

    @staticmethod
    def file_format_converter(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Converts a multiple sequence file from one format to another.

                Acceptable file formats include FASTA, Clustal, MAF, Mauve,
                Phylip, Phylip-sequential, Phylip-relaxed, and Stockholm.
                Input and output file formats are specified with the
                --input_file_format and --output_file_format arguments; input
                and output files are specified with the --input_file and
                --output_file arguments.

                Aliases:
                  file_format_converter, format_converter, ffc
                Command line interfaces: 
                  bk_file_format_converter, bk_format_converter, bk_ffc
                  

                Usage:
                biokit file_format_converter -i/--input_file <input_file>
                -iff/--input_file_format <input_file_format> 
                -o/--output_file <output_file>
                -off/--output_file_format <output_file_format>

                Options
                =====================================================
                -i/--input_file             input file name 

                -iff/--input_file_format    input file format

                -o/--output_file            output file name

                -off/--output_file_format   output file format

                Input and output file formats are specified using one of
                the following strings: fasta, clustal, maf, mauve, phylip,
                phylip_sequential, phylip_relaxed, & stockholm.
                """
            ),
        )
        parser.add_argument("-i", "--input_file", type=str, help=SUPPRESS)
        parser.add_argument("-off", "--output_file_format", type=str, help=SUPPRESS)
        parser.add_argument("-iff", "--input_file_format", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output_file", type=str, help=SUPPRESS)

        args = parser.parse_args(argv)
        FileFormatConverter(args).run()

    @staticmethod
    def gc_content(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculate GC content of a fasta file.
                
                Aliases:
                  gc_content, gc
                Command line interfaces: 
                  bk_gc_content, bk_gc

                Usage:
                biokit gc_content <fasta> [-v/--verbose]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
            
                -v, --verbose               optional argument to print
                                            the GC content of each fasta
                                            entry
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        GCContent(args).run()
    
    @staticmethod
    def gc_content_third_position(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculate GC content of the third codon position.
                The input must be the coding sequence of a gene or
                genes. All genes are assumed to have sequence lengths
                divisible by three.
                
                Aliases:
                  gc_content_third_position, gc3
                Command line interfaces: 
                  bk_gc_content_third_position, bk_gc3

                Usage:
                biokit gc_content_third_position <fasta> [-v/--verbose]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
            
                -v, --verbose               optional argument to print
                                            the GC content of each fasta
                                            entry
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        GCContentThirdPosition(args).run()

    @staticmethod
    def genome_assembly_metrics(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculate L50, L90, N50, N90, GC content, total
                length, number of scaffolds, number and sum length
                of large scaffolds, frequency of A, T, C, and G.
                
                Aliases:
                  genome_assembly_metrics, assembly_metrics
                Command line interfaces: 
                  bk_genome_assembly_metrics, bk_assembly_metrics

                Usage:
                biokit genome_assembly_metrics <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 

                -t/--threshold              threshold for what is considered
                                            a large scaffold. Only scaffolds
                                            with a length greater than this
                                            value will be counted.
                                            Default: 500
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-t","--threshold", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        GenomeAssemblyMetrics(args).run()

    @staticmethod
    def l50(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculates L50 for a genome assembly.
                
                Aliases:
                  l50
                Command line interfaces: 
                  bk_l50

                Usage:
                biokit l50 <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        L50(args).run()

    @staticmethod
    def l90(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculates L90 for a genome assembly.
                
                Aliases:
                  l90
                Command line interfaces: 
                  bk_l90

                Usage:
                biokit l90 <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        L90(args).run()

    @staticmethod
    def longest_scaffold(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate sequence length of each FASTA entry.
                
                Aliases:
                  longest_scaffold, longest_scaff, longest_contig, longest_cont
                Command line interfaces: 
                  bk_longest_scaffold, bk_longest_scaff, bk_longest_contig, bk_longest_cont

                Usage:
                biokit longest_scaffold <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        LongestScaffold(args).run()

    @staticmethod
    def n50(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculates N50 for a genome assembly.
                
                Aliases:
                  n50
                Command line interfaces: 
                  bk_n50

                Usage:
                biokit n50 <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        N50(args).run()

    @staticmethod
    def n90(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculates N90 for a genome assembly.
                
                Aliases:
                  n90
                Command line interfaces: 
                  bk_n90

                Usage:
                biokit n90 <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        N90(args).run()

    @staticmethod
    def number_of_large_scaffolds(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate number and total sequence length of
                large scaffolds. Each value is represented as
                column 1 and column 2 in the output, respectively.
                
                Aliases:
                  number_of_large_scaffolds, num_of_lrg_scaffolds,
                  number_of_large_contigs, num_of_lrg_cont
                Command line interfaces: 
                  bk_number_of_large_scaffolds, bk_num_of_lrg_scaffolds,
                  bk_number_of_large_contigs, bk_num_of_lrg_cont

                Usage:
                biokit number_of_large_scaffolds <fasta> [-t/--threshold <int>]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file
                
                -t/--threshold              threshold for what is considered
                                            a large scaffold. Only scaffolds
                                            with a length greater than this
                                            value will be counted.
                                            Default: 500
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-t","--threshold", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        NumberOfLargeScaffolds(args).run()

    @staticmethod
    def number_of_scaffolds(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate sequence length of each FASTA entry.
                
                Aliases:
                  number_of_scaffolds, num_of_scaffolds, number_of_contigs, num_of_cont
                Command line interfaces: 
                  bk_number_of_scaffolds, bk_num_of_scaffolds, bk_number_of_contigs, bk_num_of_cont

                Usage:
                biokit number_of_scaffolds <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        NumberOfScaffolds(args).run()

    @staticmethod
    def position_specific_score_matrix(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Generates a position specific score matrix for an alignment.
                
                Aliases:
                  position_specific_score_matrix, pssm
                Command line interfaces: 
                  bk_position_specific_score_matrix, bk_pssm

                Usage:
                biokit position_specific_score_matrix <fasta> [-r/--reverse]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -ac/--ambiguous_character   the ambiguity character to
                                            use. Default is 'N'
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-ac","--ambiguous_character", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        PositionSpecificScoreMatrix(args).run()

    @staticmethod
    def rename_fasta_entries(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Renames fasta entries.

                Renaming fasta entries will follow the scheme of a tab-delimited
                file wherein the first column is the current fasta entry name and
                the second column is the new fasta entry name in the resulting 
                output alignment. 

                Aliases:
                  rename_fasta_entries, rename_fasta
                Command line interfaces: 
                  bk_rename_fasta_entries, bk_rename_fasta
                
                Usage:
                biokit rename_fasta_entries <fasta> -i/--idmap <idmap>
                    [-o/--output <output_file>]
                
                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file
                -i/--idmap                  identifier map of current FASTA
                                            names (col1) and desired FASTA
                                            names (col2)
                -o/--output                 optional argument to write
                                            the renamed fasta file to.
                                            Default output has the same 
                                            name as the input file with
                                            the suffix ".renamed.fa" added
                                            to it.
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-i","--idmap", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        RenameFastaEntries(args).run()
    
    @staticmethod
    def reorder_by_sequence_length(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Reorder FASTA file entries from the longest entry
                to the shortest entry. 

                Aliases:
                  reorder_by_sequence_length, reorder_by_seq_len
                Command line interfaces: 
                  bk_reorder_by_sequence_length, bk_reorder_by_seq_len
                
                Usage:
                biokit reorder_by_sequence_length <fasta> [-o/--output <output_file>]
                
                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -o/--output                 optional argument to write
                                            the reordered fasta file to.
                                            Default output has the same 
                                            name as the input file with
                                            the suffix ".reordered.fa" added
                                            to it.
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        ReorderBySequenceLength(args).run()

    @staticmethod
    def sequence_length(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate sequence length of each FASTA entry.
                
                Aliases:
                  sequence_length, seq_len
                Command line interfaces: 
                  bk_sequence_length, bk_seq_len

                Usage:
                biokit sequence_length <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        SequenceLength(args).run()

    @staticmethod
    def sequence_complement(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Generates the sequence complement for all entries
                in a multi-FASTA file. To generate a reverse sequence
                complement, add the -r/--reverse argument.
                
                Aliases:
                  sequence_complement, seq_comp
                Command line interfaces: 
                  bk_sequence_complement, bk_seq_comp

                Usage:
                biokit sequence_complement <fasta> [-r/--reverse]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -r/--reverse                if used, the reverse complement
                                            sequence will be generated
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-r", "--reverse", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        SequenceComplement(args).run()

    @staticmethod
    def sum_of_scaffold_lengths(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Determine the sum of scaffold lengths. The
                intended use of this function is to determine
                the length of a genome assembly, but can also be
                used, for example, to determine the sum length
                of all coding sequences.
                
                Aliases:
                  sum_of_scaffold_lengths, sum_of_scaff_lens
                Command line interfaces: 
                  bk_sequence_complement, bk_sum_of_scaff_lens

                Usage:
                biokit sum_of_scaffold_lengths <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        SumOfScaffoldLengths(args).run()

    @staticmethod
    def translate_sequence(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Translates coding sequences to amino acid
                sequences. Sequences can be translated using
                diverse genetic codes. For codons that can
                encode two amino acids (e.g., TAG encodes
                Glu or STOP in the Blastocrithidia Nuclear Code),
                the standard genetic code is used.

                Custom genetic codes can be used as input and should
                be formatted with the codon in first column and the 
                resulting amino acid in the second column.
                
                Aliases:
                  translate_sequence, translate_seq, trans_seq
                Command line interfaces: 
                  bk_translate_sequence, bk_translate_seq, bk_trans_seq

                Usage:
                biokit translate_sequence <fasta> [-tt/--translation_table <code>
                -o/--output <output_file>]


                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -tt/--translation_table     Code for the translation table
                                            to be used. Default: 1, which
                                            is the standard code.
                
                -o/--output                 optional argument to write
                                            the reordered fasta file to.
                                            Default output has the same 
                                            name as the input file with
                                            the suffix ".reordered.fa" added
                                            to it.


                Codes for which translation table to use
                =====================================================
                1. The Standard Code
                2. The Vertebrate Mitochondrial Code
                3. The Yeast Mitochondrial Code
                4. The Mold, Protozoan, and Coelenterate Mitochondrial
                   Code and the Mycoplasma/Spiroplasma Code
                5. The Invertebrate Mitochondrial Code
                6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
                9. The Echinoderm and Flatworm Mitochondrial Code
                10. The Euplotid Nuclear Code
                11. The Bacterial, Archaeal and Plant Plastid Code
                12. The Alternative Yeast Nuclear Code
                13. The Ascidian Mitochondrial Code
                14. The Alternative Flatworm Mitochondrial Code
                16. Chlorophycean Mitochondrial Code
                21. Trematode Mitochondrial Code
                22. Scenedesmus obliquus Mitochondrial Code
                23. Thraustochytrium Mitochondrial Code
                24. Rhabdopleuridae Mitochondrial Code
                25. Candidate Division SR1 and Gracilibacteria Code
                26. Pachysolen tannophilus Nuclear Code
                27. Karyorelict Nuclear Code
                28. Condylostoma Nuclear Code
                29. Mesodinium Nuclear Code
                30. Peritrich Nuclear Code
                31. Blastocrithidia Nuclear Code
                33. Cephalodiscidae Mitochondrial UAA-Tyr Code
                50. CUG-Ala Code

                More information about genetic codes can be obtained from NCBI:
                https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes.
                The only codon table not described by NCBI is 50, CUG-Ala wherein CUG encodes
                for alanine.
                """
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-tt", "--translation_table", type=str, required=False, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        TranslateSequence(args).run()

    @staticmethod
    def trim_se_fastq(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Quality trim single-end FastQ data.
                
                Aliases:
                  sum_of_scaffold_lengths, sum_of_scaff_lens
                Command line interfaces: 
                  bk_sequence_complement, bk_sum_of_scaff_lens

                Usage:
                biokit sum_of_scaffold_lengths <fasta> [-m/--minimum N]

                Options
                =====================================================
                <fastq>                     first argument after 
                                            function name should be
                                            a fastq file

                -m/--minimum                minimum quality of read 
                                            to be kept (Default: 20)
                
                -l/--length                 minimum length of read 
                                            to be kept (Default: 20)

                -o/--output_file            output file name
                """
            ),
        )

        parser.add_argument("fastq", type=str, help=SUPPRESS)
        parser.add_argument("-m", "--minimum", type=str, required=False, help=SUPPRESS)
        parser.add_argument("-l", "--length", type=str, required=False, help=SUPPRESS)
        parser.add_argument("-o", "--output_file", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        TrimSEFastQ(args).run()

    ## Tree functions
    @staticmethod
    def tip_labels(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Prints the tip labels (or names) a phylogeny.

                Aliases:
                  tip_labels, tree_labels; labels; tl
                Command line interfaces: 
                  bk_tip_labels, bk_tree_labels; bk_labels; bk_tl

                Usage:
                biokit tip_labels <tree>

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        TipLabels(args).run()

def main(argv=None):
    BioKIT()

# Alignment-based functions
def faidx(argv=None):
    Biokit.faidx(sys.argv[1:])


# Tree-based functions
def tip_labels(argv=None):
    Biokit.tip_labels(sys.argv[1:])

