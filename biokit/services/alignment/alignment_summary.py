import collections
import textwrap
from typing import Any

from .base import Alignment
from ...helpers.files import read_alignment_alignio
from ...helpers.calculate_character_frequencies import calculate_character_frequencies


class AlignmentSummary(Alignment):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_format = self.normalize_output_format(self.output_format)
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
        pis, vs, cs, _ = self.determine_pis_vs_cs(alignment, alignment_length)

        # print out results
        if output_format == "tsv":
            self.print_alignment_summary(
                number_of_taxa,
                alignment_length,
                pis,
                vs,
                cs,
                char_freq,
            )
            return

        out = {
            "number_of_taxa": number_of_taxa,
            "alignment_length": alignment_length,
            "parsimony_informative_sites": pis,
            "variable_sites": vs,
            "constant_sites": cs,
            "character_frequencies": {
                key: value for key, value in sorted(char_freq.items(), key=lambda item: item[0])
            },
        }
        print(self.format_object(out, output_format))

    def combine_gaps(self, char_freq: dict[str, int]) -> dict[str, int]:
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
        char_freq: dict[str, int],
    ) -> None:
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

    def process_args(self, args: Any) -> dict[str, str | None]:
        return dict(fasta=args.fasta, output_format=getattr(args, "format", None))
