#!/usr/bin/env python

import logging
import sys
import textwrap
from .version import __version__

from argparse import (
    ArgumentParser,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .services.alignment import (
    AlignmentSummary,
    ConsensusSequence,
    PositionSpecificScoreMatrix,
)

from .services.coding_sequences import (
    GCContentFirstPosition,
    GCContentSecondPosition,
    GCContentThirdPosition,
    RelativeSynonymousCodonUsage,
    TranslateSequence,
)

from .services.fastq import (
    FastQReadLengths,
    SubsetPEFastQReads,
    SubsetSEFastQReads,
    TrimPEFastQ,
    TrimSEFastQ,
)

from .services.genome import (
    GCContent,
    GenomeAssemblyMetrics,
    L50,
    L90,
    LongestScaffold,
    N50,
    N90,
    NumberOfLargeScaffolds,
    NumberOfScaffolds,
    SumOfScaffoldLengths,
)

from .services.text import (
    CharacterFrequency,
    Faidx,
    FileFormatConverter,
    MultipleLineToSingleLineFasta,
    RenameFastaEntries,
    ReorderBySequenceLength,
    SequenceComplement,
    SequenceLength,
    SingleLineToMultipleLineFasta,
)

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

help_header = fr"""
                 ____  _       _  _______ _______
                |  _ \(_)     | |/ /_   _|__   __|
                | |_) |_  ___ | ' /  | |    | |   
                |  _ <| |/ _ \|  <   | |    | |   
                | |_) | | (_) | . \ _| |_   | |   
                |____/|_|\___/|_|\_\_____|  |_|   
                            
                Version: {__version__}
                Citation: Steenwyk et al. 2021, CITATION INFORMATION
                
"""  # noqa

translation_table_codes = f"""
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
"""  # noqa


class Biokit(object):
    help_header = fr"""
                 ____  _       _  _______ _______ 
                |  _ \(_)     | |/ /_   _|__   __|
                | |_) |_  ___ | ' /  | |    | |   
                |  _ <| |/ _ \|  <   | |    | |   
                | |_) | | (_) | . \ _| |_   | |   
                |____/|_|\___/|_|\_\_____|  |_|   
                            
                Version: {__version__}
                Citation: Steenwyk et al. 2021, CITATION INFORMATION

    """  # noqa

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

                Commands for alignments
                =======================
                alignment_summary (alias: aln_summary)
                    - calculate summary statistics for an alignment
                consensus_sequence (alias: con_seq)
                    - create a consensus sequence from an alignment
                position_specific_score_matrix (alias: pssm)
                    - create a position specific score matrix for an alignment

                Commands for coding sequences
                =============================
                gc_content_first_position (alias: gc1)
                    - calculate the GC content of the first position
                      among coding sequences
                gc_content_second_position (alias: gc2)
                    - calculate the GC content of the second position
                      among coding sequences
                gc_content_third_position (alias: gc3)
                    - calculate the GC content of the third position
                      among coding sequences
                relative_synonymous_codon_usage (alias: rscu)
                    - calculate relative synonymous codon usage
                      to evaluate potential codon usage biases
                translate_sequence (alias: translate_seq, trans_seq) 
                    - translate coding sequences to amino acids

                Commands for fastq files
                ========================
                fastq_read_lengths (alias: fastq_read_lens)
                    - determine the lengths of fastq reads
                subset_pe_fastq_reads (alias: subset_pe_fastq)
                    - subset paired-end fastq reads and
                      maintain pairing information
                subset_se_fastq_reads (alias: subset_se_fastq)
                    - subset single-end fastq reads
                trim_pe_fastq
                    - quality trim paired-end fastq reads
                      and maintain pairing information
                trim_se_fastq
                    - quality trim single-end fastq reads

                Commands for genomes
                ====================
                gc_content (alias: gc)
                    - calculate the GC content of a FASTA file
                genome_assembly_metrics (alias: assembly_metrics)
                    - calculate various genome assembly metrics
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
                sum_of_scaffold_lengths (alias: sum_of_contig_lengths)
                    - calculate sum of scaffold/contig lengths

                Commands for sequence files
                ===========================
                character_frequency (alias: char_freq)
                    - determine the frequency of all observed characters
                faidx (alias: get_entry; ge)
                    - extract query fasta entry from multi-fasta file
                file_format_converter (alias: format_converter; ffc)
                    - convert a multiple sequence file from one format
                      to another
                multiple_line_to_single_line_fasta (alias: ml2sl)
                    - reformats sequences that occur on multiple
                      lines to be represented in a single line
                rename_fasta_entries (alias: rename_fasta)
                    - rename entries in a FASTA file
                reorder_by_sequence_length (alias: reorder_by_seq_len)
                    - reorder sequences from longest to shortest in a FASTA file
                sequence_complement (alias: seq_comp)
                    - generate the complementary sequence for an alignment 
                sequence_length (alias: seq_len)
                    - calculate the length of each FASTA entry
                single_line_to_multiple_line_fasta (alias: sl2ml)
                    - reformats sequences so that there are 60
                      characters per sequence line
                """  # noqa
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
                self.run_alias(args.command, sys.argv[2:], parser)
        except NameError as e:
            print(e)
            sys.exit()

    # aliases # noqa
    def run_alias(self, command, argv, parser):
        # version
        if command in ["v"]:
            return self.version()
        # aliases for alignments
        elif command in ["aln_summary"]:
            return self.alignment_summary(argv)
        elif command in ["con_seq"]:
            return self.consensus_sequence(argv)
        elif command in ["pssm"]:
            return self.position_specific_score_matrix(argv)
        # aliases for coding sequences
        elif command in ["gc1"]:
            return self.gc_content_first_position(argv)
        elif command in ["gc2"]:
            return self.gc_content_second_position(argv)
        elif command in ["gc3"]:
            return self.gc_content_third_position(argv)
        elif command in ["rscu"]:
            return self.relative_synonymous_codon_usage(argv)
        elif command in ["translate_seq", "trans_seq"]:
            return self.translate_sequence(argv)
        # aliases for fastq files
        elif command in ["fastq_read_lens"]:
            return self.fastq_read_lengths(argv)
        elif command in ["subset_pe_fastq"]:
            return self.subset_pe_fastq_reads(argv)
        elif command in ["subset_se_fastq"]:
            return self.subset_se_fastq_reads(argv)
        # aliases for genomes
        elif command in ["gc"]:
            return self.gc_content(argv)
        elif command in ["assembly_metrics"]:
            return self.genome_assembly_metrics(argv)
        elif command in ["longest_scaff", "longest_contig", "longest_cont"]:
            return self.longest_scaffold(argv)
        elif command in [
            "num_of_lrg_scaffolds",
            "number_of_large_contigs",
            "num_of_lrg_cont",
        ]:
            return self.number_of_large_scaffolds(argv)
        elif command in ["num_of_scaffolds", "number_of_contigs", "num_of_cont"]:
            return self.number_of_scaffolds(argv)
        elif command in ["sum_of_contig_lengths"]:
            return self.sum_of_scaffold_lengths(argv)
        # alias for sequence files
        elif command in ["char_freq"]:
            return self.character_frequency(argv)
        elif command in ["get_entry", "ge"]:
            return self.faidx(argv)
        elif command in ["format_converter", "ffc"]:
            return self.file_format_converter(argv)
        elif command in ["ml2sl"]:
            return self.multiple_line_to_single_line_fasta(argv)
        elif command in ["rename_fasta"]:
            return self.rename_fasta_entries(argv)
        elif command in ["reorder_by_seq_len"]:
            return self.reorder_by_sequence_length(argv)
        elif command in ["seq_comp"]:
            return self.sequence_complement(argv)
        elif command in ["seq_len"]:
            return self.sequence_length(argv)
        elif command in ["sl2ml"]:
            return self.single_line_to_multiple_line_fasta(argv)
        else:
            print(
                "Invalid command option. See help for a complete list of commands and aliases."
            )
            parser.print_help()
            sys.exit(1)

    # print version
    def version(self):
        print(
            textwrap.dedent(
                f"""\
            {self.help_header}
            """
            )
        )

    # alignment functions
    @staticmethod
    def alignment_summary(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Summary statistics for an alignment. Reported
                statistics include alignment length, number of taxa,
                number of parsimony sites, number of variable sites,
                number of constant sites, frequency of each character
                (including gaps, which are considered to be '-' or '?'). 
                
                Aliases:
                  alignment_summary, aln_summary
                Command line interfaces: 
                  bk_alignment_summary, bk_aln_summary

                Usage:
                biokit alignment_summary <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        AlignmentSummary(args).run()

    @staticmethod
    def consensus_sequence(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Generates a consequence from a multiple sequence alignment
                file in FASTA format.
                
                Aliases:
                  consensus_sequence, con_seq
                Command line interfaces: 
                  bk_consensus_sequence, bk_con_seq

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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-t", "--threshold", type=str, help=SUPPRESS)
        parser.add_argument("-ac", "--ambiguous_character", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        ConsensusSequence(args).run()

    @staticmethod
    def position_specific_score_matrix(argv):
        parser = ArgumentParser(
            add_help=True,
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
                biokit position_specific_score_matrix <fasta> 
                [-ac/--ambiguous_character <ambiguous character>]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -ac/--ambiguous_character   the ambiguity character to
                                            use. Default is 'N'
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-ac", "--ambiguous_character", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        PositionSpecificScoreMatrix(args).run()

    # coding sequence functions
    @staticmethod
    def gc_content_first_position(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculate GC content of the first codon position.
                The input must be the coding sequence of a gene or
                genes. All genes are assumed to have sequence lengths
                divisible by three.
                
                Aliases:
                  gc_content_first_position, gc1
                Command line interfaces: 
                  bk_gc_content_first_position, bk_gc1

                Usage:
                biokit gc_content_first_position <fasta> [-v/--verbose]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
            
                -v, --verbose               optional argument to print
                                            the GC content of each fasta
                                            entry
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        GCContentFirstPosition(args).run()

    @staticmethod
    def gc_content_second_position(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculate GC content of the second codon position.
                The input must be the coding sequence of a gene or
                genes. All genes are assumed to have sequence lengths
                divisible by three.
                
                Aliases:
                  gc_content_second_position, gc2
                Command line interfaces: 
                  bk_gc_content_second_position, bk_gc2

                Usage:
                biokit gc_content_second_position <fasta> [-v/--verbose]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
            
                -v, --verbose               optional argument to print
                                            the GC content of each fasta
                                            entry
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        GCContentSecondPosition(args).run()

    @staticmethod
    def gc_content_third_position(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        GCContentThirdPosition(args).run()

    @staticmethod
    def relative_synonymous_codon_usage(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Calculate relative synonymous codon usage.

                Relative synonymous codon usage is the ratio
                of the observed frequency of codons over the
                expected frequency given that all the synonymous
                codons for the same amino acids are used equally.

                Aliases:
                  relative_synonymous_codon_usage, rscu
                Command line interfaces: 
                  bk_relative_synonymous_codon_usage, bk_rscu
                
                Usage:
                biokit relative_synonymous_codon_usage <fasta> 
                
                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -tt/--translation_table     Code for the translation table
                                            to be used. Default: 1, which
                                            is the standard code.


                {translation_table_codes}
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument(
            "-tt", "--translation_table", type=str, required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        RelativeSynonymousCodonUsage(args).run()

    @staticmethod
    def translate_sequence(argv):
        parser = ArgumentParser(
            add_help=True,
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


                {translation_table_codes}
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument(
            "-tt", "--translation_table", type=str, required=False, help=SUPPRESS
        )
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        TranslateSequence(args).run()

    # fastq file functions
    @staticmethod
    def fastq_read_lengths(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Quality trim single-end FastQ data.
                
                Aliases:
                  fastq_read_lengths, fastq_read_lens
                Command line interfaces: 
                  bk_fastq_read_lengths, bk_fastq_read_lens

                Usage:
                biokit fastq_read_lengths <fasta> [-v/--verbose]

                Options
                =====================================================
                <fastq>                     first argument after 
                                            function name should be
                                            a fastq file

                -v/--verbose                print length of each fastq
                                            read
                """  # noqa
            ),
        )

        parser.add_argument("fastq", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        FastQReadLengths(args).run()

    @staticmethod
    def subset_pe_fastq_reads(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Subset paired-end FastQ data.
                
                Aliases:
                  subset_pe_fastq_reads, subset_pe_fastq
                Command line interfaces: 
                  bk_subset_pe_fastq_reads, bk_subset_pe_fastq

                Usage:
                biokit subset_pe_fastq_reads <fasta>

                Options
                =====================================================
                <fastq1>                    first argument after 
                                            function name should be
                                            a fastq file

                <fastq2>                    second argument after 
                                            function name should be
                                            a fastq file

                -p/--percent                percentage of reads to
                                            maintain in subsetted data.
                                            Default: 10

                -s/--seed                   seed for random sampling.
                                            Default: date and time

                -o/--output_file            output file name
                """  # noqa
            ),
        )

        parser.add_argument("fastq1", type=str, help=SUPPRESS)
        parser.add_argument("fastq2", type=str, help=SUPPRESS)
        parser.add_argument("-p", "--percent", type=str, required=False, help=SUPPRESS)
        parser.add_argument("-s", "--seed", type=str, required=False, help=SUPPRESS)
        parser.add_argument(
            "-o", "--output_file", type=str, required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        SubsetPEFastQReads(args).run()

    @staticmethod
    def subset_se_fastq_reads(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Subset single-end FastQ data.
                
                Aliases:
                  subset_se_fastq_reads, subset_se_fastq
                Command line interfaces: 
                  bk_subset_se_fastq_reads, bk_subset_se_fastq

                Usage:
                biokit subset_se_fastq_reads <fasta>

                Options
                =====================================================
                <fastq>                     first argument after 
                                            function name should be
                                            a fastq file

                -p/--percent                percentage of reads to
                                            maintain in subsetted data.
                                            Default: 10

                -s/--seed                   seed for random sampling.
                                            Default: date and time

                -o/--output_file            output file name
                """  # noqa
            ),
        )

        parser.add_argument("fastq", type=str, help=SUPPRESS)
        parser.add_argument("-p", "--percent", type=str, required=False, help=SUPPRESS)
        parser.add_argument("-s", "--seed", type=str, required=False, help=SUPPRESS)
        parser.add_argument(
            "-o", "--output_file", type=str, required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        SubsetSEFastQReads(args).run()

    @staticmethod
    def trim_pe_fastq(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Quality trim paired-end FastQ data.
                
                Aliases:
                  trim_pe_fastq
                Command line interfaces: 
                  bk_trim_pe_fastq

                Usage:
                biokit trim_pe_fastq <fastq1> <fastq2> [-m/--minimum 20
                    -l/--length 20]

                Options
                =====================================================
                <fastq1>                    first argument after 
                                            function name should be
                                            a fastq file

                <fastq2>                    second argument after 
                                            function name should be
                                            a fastq file

                -m/--minimum                minimum quality of read 
                                            to be kept (Default: 20)
                
                -l/--length                 minimum length of read 
                                            to be kept (Default: 20)
                """  # noqa
            ),
        )

        parser.add_argument("fastq1", type=str, help=SUPPRESS)
        parser.add_argument("fastq2", type=str, help=SUPPRESS)
        parser.add_argument("-m", "--minimum", type=str, required=False, help=SUPPRESS)
        parser.add_argument("-l", "--length", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        TrimPEFastQ(args).run()

    @staticmethod
    def trim_se_fastq(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fastq", type=str, help=SUPPRESS)
        parser.add_argument("-m", "--minimum", type=str, required=False, help=SUPPRESS)
        parser.add_argument("-l", "--length", type=str, required=False, help=SUPPRESS)
        parser.add_argument(
            "-o", "--output_file", type=str, required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        TrimSEFastQ(args).run()

    # genome functions
    @staticmethod
    def gc_content(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        GCContent(args).run()

    @staticmethod
    def genome_assembly_metrics(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-t", "--threshold", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        GenomeAssemblyMetrics(args).run()

    @staticmethod
    def l50(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        L50(args).run()

    @staticmethod
    def l90(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        L90(args).run()

    @staticmethod
    def longest_scaffold(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        LongestScaffold(args).run()

    @staticmethod
    def n50(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        N50(args).run()

    @staticmethod
    def n90(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        N90(args).run()

    @staticmethod
    def number_of_large_scaffolds(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-t", "--threshold", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        NumberOfLargeScaffolds(args).run()

    @staticmethod
    def number_of_scaffolds(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        NumberOfScaffolds(args).run()

    @staticmethod
    def sum_of_scaffold_lengths(argv):
        parser = ArgumentParser(
            add_help=True,
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
                  sum_of_scaffold_lengths, sum_of_contig_lengths
                Command line interfaces: 
                  bk_sequence_complement, bk_sum_of_contig_lengths

                Usage:
                biokit sum_of_scaffold_lengths <fasta>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        SumOfScaffoldLengths(args).run()

    # text functions
    @staticmethod
    def character_frequency(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        CharacterFrequency(args).run()

    @staticmethod
    def faidx(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-e", "--entry", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        Faidx(args).run()

    @staticmethod
    def file_format_converter(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("-i", "--input_file", type=str, help=SUPPRESS)
        parser.add_argument("-off", "--output_file_format", type=str, help=SUPPRESS)
        parser.add_argument("-iff", "--input_file_format", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output_file", type=str, help=SUPPRESS)

        args = parser.parse_args(argv)
        FileFormatConverter(args).run()

    @staticmethod
    def multiple_line_to_single_line_fasta(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Converts FASTA files with multiple lines
                per sequence to a FASTA file with the sequence
                on one line.

                Aliases:
                  multiple_line_to_single_line_fasta, ml2sl
                Command line interfaces: 
                  bk_multiple_line_to_single_line_fasta, bk_ml2sl
                
                Usage:
                biokit multiple_line_to_single_line_fasta <fasta> 
                    [-o/--output <output_file>]
                
                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        MultipleLineToSingleLineFasta(args).run()

    @staticmethod
    def rename_fasta_entries(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-i", "--idmap", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        RenameFastaEntries(args).run()

    @staticmethod
    def reorder_by_sequence_length(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        ReorderBySequenceLength(args).run()

    @staticmethod
    def sequence_complement(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument(
            "-r", "--reverse", action="store_true", required=False, help=SUPPRESS
        )
        args = parser.parse_args(argv)
        SequenceComplement(args).run()

    @staticmethod
    def sequence_length(argv):
        parser = ArgumentParser(
            add_help=True,
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
                """  # noqa
            ),
        )

        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        SequenceLength(args).run()

    @staticmethod
    def single_line_to_multiple_line_fasta(argv):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Converts FASTA files with single lines per
                sequence to a FASTA file with the sequence
                on multiple lines. Each line with have 60 
                characters following standard NCBI format.

                Aliases:
                  single_line_to_multiple_line_fasta, sl2ml
                Command line interfaces: 
                  bk_single_line_to_multiple_line_fasta, bk_sl2ml
                
                Usage:
                biokit single_line_to_multiple_line_fasta <fasta>
                
                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file
                """  # noqa
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        SingleLineToMultipleLineFasta(args).run()


def main(argv=None):
    Biokit()


# Alignment-based functions
def faidx(argv=None):
    Biokit.faidx(sys.argv[1:])
