import sys

from Bio import SeqIO

from .base import Text

class RelativeSynonymousCodonUsage(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        translation_table = self.read_translation_table()

        # get codon_table
        codon_table = dict()
        for k, v in translation_table.items():
            if v in codon_table.keys():
                codon_table[v].append(k)
            else:
                codon_table[v] = [k]

        # get counts of codons
        codon_counts = dict()
        records = SeqIO.parse(self.fasta, "fasta")
        for seq_record in records:
            if len(seq_record._seq) % 3 == 0:
                for position in range(0, len(seq_record._seq), 3):
                    codon = seq_record._seq[position:position+3]._data.upper().replace("T", "U")
                    if codon in codon_counts.keys():
                        codon_counts[codon] += 1
                    else:
                        codon_counts[codon] = 1
        
        # calculate rscu
        rscu = dict()
        for aa, codons in codon_table.items():
            observed_sum = 0
            for codon in codons:
                observed_sum += codon_counts[codon]
            for codon in codons:
                rscu[codon] = round(codon_counts[codon]/(observed_sum/len(codons)), 4)
        
        # reverse sort according to rscu values
        rscu = dict(sorted(rscu.items(), key=lambda item: item[1], reverse=True))
        
        for codon, rscu in rscu.items():
            print(f"{codon}\t{rscu}")
    
    def process_args(self, args):
        return dict(
            fasta=args.fasta,
            translation_table=args.translation_table
        )