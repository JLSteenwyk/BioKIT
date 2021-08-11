import sys

from Bio import SeqIO

from .base import Text

class TranslateSequence(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get the translation table as a dictionary
        translation_table = self.read_translation_table(self.translation_table)

        records = SeqIO.parse(self.fasta, "fasta")
        
        for seq_record in records:
            amino_acids = []
            if len(seq_record._seq) % 3 == 0:
                for position in range(0, len(seq_record._seq), 3):
                    codon = seq_record._seq[position:position+3]._data.upper().replace("T", "U")
                    amino_acids.append(translation_table[codon])
            
            print(f">{seq_record.id}")
            print('\n'.join(''.join(amino_acids)[i:i+80] for i in range(0, len(amino_acids), 80)))
    
    def read_translation_table(self, translation_table: str):
        translation_table = dict()
        if self.translation_table is None:
            pathing = '/Users/jlsteenwyk/Desktop/GITHUB/BioKIT/biokit/tables/standard_genetic_code.txt'
            with open(pathing) as code:
                for line in code:
                    line=line.split()
                    translation_table[line[0]] = line[1]
        elif self.translation_table in [
            '1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14',
            '16', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30',
            '31', '33', '50'
        ]:
            print("eh")
        # case handling for a custom translation table
        else:
            with open(self.translation_table) as code:
                for line in code:
                    line=line.split()
                    translation_table[line[0]] = line[1]
        
        return translation_table


    def process_args(self, args):
        return dict(fasta=args.fasta, translation_table=args.translation_table)