import collections
import sys
import textwrap

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

        if "?" in char_freq.keys() and "-" in char_freq.keys():
            char_freq["-"] += char_freq["?"]
            del char_freq["?"]

        # number of taxa
        number_of_taxa = len(alignment)

        # alignment length
        alignment_length = alignment.get_alignment_length()

        # number of parsimony informative sites
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

        print(
            textwrap.dedent(
                f"\
            General Characteristics\n\
            =======================\n\
            {number_of_taxa}\tNumber of taxa\n\
            {alignment_length}\tAlignment length\n\
            {parsimony_informative_sites}\tParsimony informative sites\n\
            {variable_sites}\tVariable sites\n\
            {constant_sites}\tConstant sites"
            )
        )
        char_freq = collections.OrderedDict(sorted(char_freq.items(), reverse=True))
        print("\nCharacter Frequencies")
        print("=====================")
        for k, v in char_freq.items():
            print(f"{k}\t{v}")

        # Summary statistics for an alignment. Reported
        # statistics include alignment length, number of taxa,
        # number of parsimony sites, number of variable sites,
        # number of constant sites, frequency of each character
        # (including gaps, which are considered to be '-' or '?').

    def process_args(self, args):
        return dict(fasta=args.fasta)
