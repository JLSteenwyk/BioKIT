from os import path
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .base import Text

here = path.dirname(__file__)

class TranslateSequence(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get the translation table as a dictionary
        translation_table = self.read_translation_table(self.translation_table)

        records = SeqIO.parse(self.fasta, "fasta")
        
        with open(self.output_file_path, 'w') as output_file_path:
            for seq_record in records:
                amino_acids = []
                if len(seq_record._seq) % 3 == 0:
                    for position in range(0, len(seq_record._seq), 3):
                        codon = seq_record._seq[position:position+3]._data.upper().replace("T", "U")
                        amino_acids.append(translation_table[codon])

                translated_seq_record = SeqRecord(
                    Seq(''.join(amino_acids)),
                    id=seq_record.id,
                    name='',
                    description='',
                )

                SeqIO.write(translated_seq_record, output_file_path, "fasta")
    
    def read_translation_table(self, translation_table: str):
        trans_table = dict()

        if translation_table is None or translation_table == '1':
            pathing = path.join(here, "../../tables/standard_genetic_code.txt")
        elif translation_table == '2':
            pathing = path.join(here, "../../tables/vertebrate_mitochondrial_code.txt")
        elif translation_table == '3':
            pathing = path.join(here, "../../tables/yeast_mitochondrial_code.txt")
        elif translation_table == '4':
            pathing = path.join(here, "../../tables/mold_protozoan_and_coelenterate_mitochondrial_code_and_the_mycoplasma_spiroplasma.txt")
        elif translation_table == '5':
            pathing = path.join(here, "../../tables/invertebrate_mitochondrial_code.txt")
        elif translation_table == '6':
            pathing = path.join(here, "../../tables/ciliate_dasycladacean_and_hexamita_nuclear_code.txt")
        elif translation_table == '9':
            pathing = path.join(here, "../../tables/echinoderm_and_flatworm_mitochondrial_code.txt")
        elif translation_table == '10':
            pathing = path.join(here, "../../tables/euplotid_nuclear_code.txt")
        elif translation_table == '11':
            pathing = path.join(here, "../../tables/bacterial_archaeal_and_plant_plastid_code.txt")
        elif translation_table == '12':
            pathing = path.join(here, "../../tables/alternative_yeast_nuclear_code.txt")
        elif translation_table == '13':
            pathing = path.join(here, "../../tables/ascidian_mitochondrial_code.txt")
        elif translation_table == '14':
            pathing = path.join(here, "../../tables/alternative_flatworm_mitochondrial_code.txt")
        elif translation_table == '16':
            pathing = path.join(here, "../../tables/chlorophycean_mitochondrial_code.txt")
        elif translation_table == '21':
            pathing = path.join(here, "../../tables/trematode_mitochondrial_code.txt")
        elif translation_table == '22':
            pathing = path.join(here, "../../tables/scenedesmus_obliquus_mitochondrial_code.txt")
        elif translation_table == '23':
            pathing = path.join(here, "../../tables/thraustochytrium_mitochondrial_code.txt")
        elif translation_table == '24':
            pathing = path.join(here, "../../tables/rhabdopleuridae_mitochondrial_code.txt")
        elif translation_table == '25':
            pathing = path.join(here, "../../tables/candidate_division_sr1_and_gracilibacteria_code.txt")
        elif translation_table == '26':
            pathing = path.join(here, "../../tables/pachysolen_tannophilus_nuclear_code.txt")
        elif translation_table == '27':
            pathing = path.join(here, "../../tables/karyorelict_nuclear_code.txt")
        elif translation_table == '28':
            pathing = path.join(here, "../../tables/condylostoma_nuclear_code.txt")
        elif translation_table == '29':
            pathing = path.join(here, "../../tables/mesodinium_nuclear_code.txt")
        elif translation_table == '30':
            pathing = path.join(here, "../../tables/peritrich_nuclear_code.txt")
        elif translation_table == '31':
            pathing = path.join(here, "../../tables/blastocrithidia_nuclear_code.txt")
        elif translation_table == '33':
            pathing = path.join(here, "../../tables/cephalodiscidae_mitochondrial_UAA_tyr_code.txt")
        elif translation_table == '50':
            pathing = path.join(here, "../../tables/CUG_ala_code.txt")
        # case handling for a custom translation table
        else:
            trans_table = self.read_lookup_table(translation_table, trans_table)
        
        return self.read_lookup_table(pathing, trans_table)

    def read_lookup_table(self, pathing: str, trans_table: dict):
        with open(pathing) as code:
            for line in code:
                line=line.split()
                trans_table[line[0]] = line[1]
        return trans_table

    def process_args(self, args):
        if args.output is None:
            output_file_path = f"{args.fasta}.translated.fa"
        else:
            output_file_path = f"{args.output}"

        return dict(
            fasta=args.fasta,
            translation_table=args.translation_table,
            output_file_path=output_file_path,
        )