import sys
import textwrap

from .base import Text

class GenomeAssemblyMetrics(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        l50 = self.calc_l50()
        l90 = self.calc_l90()
        n50 = self.calc_n50()
        n90 = self.calc_n90()
        
        char_freq = self.character_frequency()
        sum_of_chars = sum(char_freq.values())
        gc_content = (char_freq['G']+char_freq['C'])/sum_of_chars
        a_freq = char_freq['A']/sum_of_chars
        t_freq = char_freq['T']/sum_of_chars
        g_freq = char_freq['G']/sum_of_chars
        c_freq = char_freq['C']/sum_of_chars

        assembly_size = self.sum_of_scaffold_lengths()

        num_of_scaffolds = self.number_of_scaffolds()

        num_of_large_scaffolds, len_of_large_scaffolds = self.number_of_large_scaffolds()

        print(textwrap.dedent(f"\
            Assembly size\t{sum_of_chars}\n\
            L50\t{l50}\n\
            L90\t{l90}\n\
            N50\t{n50}\n\
            N90\t{n90}\n\
            GC content\t{round(gc_content, 4)}\n\
            Number of scaffolds\t{num_of_scaffolds}\n\
            Number of large scaffolds\t{num_of_large_scaffolds}\n\
            Sum length of large scaffolds\t{len_of_large_scaffolds}\n\
            Frequency of A\t{round(a_freq, 4)}\n\
            Frequency of T\t{round(t_freq, 4)}\n\
            Frequency of C\t{round(c_freq, 4)}\n\
            Frequency of G\t{round(g_freq, 4)}"
        ))
    
    def process_args(self, args):
        if args.threshold is None:
            threshold = 500
        else:
            threshold = int(args.threshold)

        return dict(fasta=args.fasta, threshold=threshold)