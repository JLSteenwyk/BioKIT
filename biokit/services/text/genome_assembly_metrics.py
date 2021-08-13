import sys

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



        print(f"{l50}\t{l90}\t{n50}\t{n90}\t{char_freq}\t{gc_content}\t{num_of_scaffolds}\t{num_of_large_scaffolds}\t{len_of_large_scaffolds}\t{a_freq}\t{t_freq}\t{c_freq}\t{g_freq}")
        # Calculate frequency of A, T, C, and G.
    
    def process_args(self, args):
        if args.threshold is None:
            threshold = 500
        else:
            threshold = int(args.threshold)

        return dict(fasta=args.fasta, threshold=threshold)