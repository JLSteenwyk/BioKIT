import collections
import textwrap
import sys

from Bio import AlignIO

from .base import Alignment
from ...helpers.calculate_character_frequencies import calculate_character_frequencies


class AlignmentSummary(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = AlignIO.read(open(self.fasta), "fasta")

        # character frequency counting
        char_freq = calculate_character_frequencies(self.fasta)
        char_freq = self.combine_gaps(char_freq)
        char_freq = collections.OrderedDict(sorted(char_freq.items(), reverse=True))

        # number of taxa
        number_of_taxa = len(alignment)

        # alignment length
        alignment_length = alignment.get_alignment_length()

        # number of parsimony informative, variable, and constant sites
        pis, vs, cs = self.determine_pis_vs_cs(alignment, alignment_length)

        # print out results
        self.print_alignment_summary(
            number_of_taxa,
            alignment_length,
            pis,
            vs,
            cs,
            char_freq,
        )

    def determine_pis_vs_cs(self, alignment, alignment_length):
        """
        determine number of parsimony informative,
        variable, and constant sites in an alignment
        """
        parsimony_informative_sites = 0
        variable_sites = 0
        constant_sites = 0
        for i in range(0, alignment_length):
            seq_at_position = ""
            seq_at_position += alignment[:, i]
            seq_at_position = seq_at_position.replace("?", "")
            seq_at_position = seq_at_position.replace("-", "")
            num_occurences = {}
            for char in set(seq_at_position):
                num_occurences[char] = seq_at_position.count(char)
            d = dict((k, v) for k, v in num_occurences.items() if v >= 2)

            # two characters that occur at least twice
            if len(d) >= 2:
                parsimony_informative_sites += 1
                variable_sites += 1
            # if one character occurs at least twice and is the only character,
            # the site is not parismony informative but it is constant
            elif len(d) == 1 and len(num_occurences) == 1:
                constant_sites += 1
            else:
                1

        return parsimony_informative_sites, variable_sites, constant_sites

    def combine_gaps(self, char_freq: dict) -> dict:
        """
        combine counts for '?' and '-' characters,
        which represents gaps
        """
        if "?" in char_freq.keys() and "-" in char_freq.keys():
            char_freq["-"] += char_freq["?"]
            del char_freq["?"]

        return char_freq

    def print_alignment_summary(
        self,
        number_of_taxa: int,
        alignment_length: int,
        parsimony_informative_sites: int,
        variable_sites: int,
        constant_sites: int,
        char_freq: dict,
    ):
        """
        print summary results
        """

        out_str = textwrap.dedent(
            f"""
            General Characteristics\n\
            =======================\n\
            {number_of_taxa}\tNumber of taxa\n\
            {alignment_length}\tAlignment length\n\
            {parsimony_informative_sites}\tParsimony informative sites\n\
            {variable_sites}\tVariable sites\n\
            {constant_sites}\tConstant sites"""
        )

        out_str += "\n\nCharacter Frequencies"
        out_str += "\n====================="
        for k, v in char_freq.items():
            out_str += f"\n{k}\t{v}"

        print(out_str)

    def process_args(self, args):
        return dict(fasta=args.fasta)
