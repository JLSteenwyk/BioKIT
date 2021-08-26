import collections
import textwrap

from .base import Alignment
from ...helpers.files import read_alignment_alignio
from ...helpers.calculate_character_frequencies import calculate_character_frequencies


class AlignmentSummary(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = read_alignment_alignio(self.fasta)

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
