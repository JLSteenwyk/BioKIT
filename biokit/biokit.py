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
    Faidx,
    FileFormatConverter,
    RenameFastaEntries,
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

                BioKIT helps process and analyze multiple sequence alignments and phylogenies.

                Generally, all functions are designed to help understand the contents of alignments
                (e.g., gc content or the number of parsimony informative sites) and the shape
                of trees (e.g., treeness, degree of violation of a molecular clock).

                Some help messages indicate that summary statistics are reported (e.g., 
                bipartition_support_stats). Summary statistics include mean, median, 25th percentile,
                75th percentile, minimum, maximum, standard deviation, and variance. These functions
                typically have a verbose option that allows users to get the underlying data
                used to calculate summary statistics. 

                Usage: phykit <command> [optional command arguments]

                Command specific help messages can be viewed by adding a 
                -h/--help argument after the command. For example, to see the
                to see the help message for the command 'treeness', execute
                "phykit treeness -h" or "phykit treeness --help".

                Lastly, each function comes with aliases to save the user some
                key strokes. For example, to get the help message for the 'treeness'
                function, you can type "phykit tness -h". All aliases are specified
                in parentheses after the long form of the function name. 

                Alignment-based commands
                ========================
                faidx (alias: get_entry; ge)
                    - extract query fasta entry from multi-fasta file                              

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
        # Alignment aliases
        elif command in ['get_entry', 'ge']:
            return self.faidx(argv)
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
                  faidx, get_entry; ge
                Command line interfaces: 
                  pk_faidx, pk_get_entry, pk_ge
                  

                Usage:
                phykit faidx <fasta> -e/--entry <fasta entry>

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

                REWRITE
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
                  pk_rename_fasta_entries, pk_rename_fasta
                Usage:
                phykit rename_fasta_entries <fasta> -i/--idmap <idmap>
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
                  pk_tip_labels, pk_tree_labels; pk_labels; pk_tl

                Usage:
                phykit tip_labels <tree>

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

