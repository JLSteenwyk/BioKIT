from os import path
import sys

from Bio import SeqIO

from .base import Text

here = path.dirname(__file__)

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
            pathing = path.join(here, "../../tables/standard_genetic_code.txt")
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
            # 1. The Standard Code
            # 2. The Vertebrate Mitochondrial Code
            # 3. The Yeast Mitochondrial Code
            # 4. The Mold, Protozoan, and Coelenterate Mitochondrial
            #     Code and the Mycoplasma/Spiroplasma Code
            # 5. The Invertebrate Mitochondrial Code
            # 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
            # 9. The Echinoderm and Flatworm Mitochondrial Code
            # 10. The Euplotid Nuclear Code
            # 11. The Bacterial, Archaeal and Plant Plastid Code
            # 12. The Alternative Yeast Nuclear Code
            # 13. The Ascidian Mitochondrial Code
            # 14. The Alternative Flatworm Mitochondrial Code
            # 16. Chlorophycean Mitochondrial Code
            # 21. Trematode Mitochondrial Code
            # 22. Scenedesmus obliquus Mitochondrial Code
            # 23. Thraustochytrium Mitochondrial Code
            # 24. Rhabdopleuridae Mitochondrial Code
            # 25. Candidate Division SR1 and Gracilibacteria Code
            # 26. Pachysolen tannophilus Nuclear Code
            # 27. Karyorelict Nuclear Code
            # 28. Condylostoma Nuclear Code
            # 29. Mesodinium Nuclear Code
            # 30. Peritrich Nuclear Code
            # 31. Blastocrithidia Nuclear Code
            # 33. Cephalodiscidae Mitochondrial UAA-Tyr Code
            # 50. CUG-Ala Code
        # case handling for a custom translation table
        else:
            with open(self.translation_table) as code:
                for line in code:
                    line=line.split()
                    translation_table[line[0]] = line[1]
        
        return translation_table


    def process_args(self, args):
        return dict(fasta=args.fasta, translation_table=args.translation_table)