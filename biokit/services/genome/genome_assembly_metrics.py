import textwrap

from .base import Genome
from ...helpers.calculate_character_frequencies import calculate_character_frequencies


class GenomeAssemblyMetrics(Genome):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        l50 = self.calc_l50()
        l90 = self.calc_l90()
        n50 = self.calc_n50()
        n90 = self.calc_n90()

        char_freq = calculate_character_frequencies(self.fasta)
        assembly_size = self.sum_of_scaffold_lengths()

        gc_content = (char_freq["G"] + char_freq["C"]) / assembly_size
        a_freq = char_freq["A"] / assembly_size
        t_freq = char_freq["T"] / assembly_size
        g_freq = char_freq["G"] / assembly_size
        c_freq = char_freq["C"] / assembly_size

        num_of_scaffolds = self.number_of_scaffolds()

        longest_scaffold = self.longest_scaffold()

        (
            num_of_large_scaffolds,
            len_of_large_scaffolds,
        ) = self.number_of_large_scaffolds()

        print(
            textwrap.dedent(
                f"\
            {assembly_size}\tAssembly size\n\
            {l50}\tL50\n\
            {l90}\tL90\n\
            {n50}\tN50\n\
            {n90}\tN90\n\
            {round(gc_content, 4)}\tGC content\n\
            {num_of_scaffolds}\tNumber of scaffolds\n\
            {num_of_large_scaffolds}\tNumber of large scaffolds\n\
            {len_of_large_scaffolds}\tSum length of large scaffolds\n\
            {longest_scaffold}\tLongest scaffold\n\
            {round(a_freq, 4)}\tFrequency of A\n\
            {round(t_freq, 4)}\tFrequency of T\n\
            {round(c_freq, 4)}\tFrequency of C\n\
            {round(g_freq, 4)}\tFrequency of G"
            )
        )

    def process_args(self, args):
        if args.threshold is None:
            threshold = 500
        else:
            threshold = int(args.threshold)

        return dict(fasta=args.fasta, threshold=threshold)
