from os import path

from Bio import SeqIO

from ..base import BaseService

here = path.dirname(__file__)


class CodingSequence(BaseService):
    def __init__(
        self,
        *args,
        fasta=None,
        verbose=None,
        translation_table=None,
        output_file_path=None,
    ):
        self.fasta = fasta
        self.verbose = verbose
        self.translation_table = translation_table
        self.output_file_path = output_file_path

    def calculate_rscu(self, translation_table) -> dict:
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
                    codon = (
                        seq_record._seq[position:position + 3]
                        ._data.upper()
                        .replace("T", "U")
                    )
                    if codon in codon_counts.keys():
                        codon_counts[codon] += 1
                    else:
                        codon_counts[codon] = 1

        # calculate rscu
        rscu = dict()
        for _, codons in codon_table.items():
            observed_sum = 0
            for codon in codons:
                try:
                    observed_sum += codon_counts[codon]
                except KeyError:
                    observed_sum = 0
            for codon in codons:
                try:
                    rscu[codon] = round(
                        codon_counts[codon] / (observed_sum / len(codons)), 4
                    )
                except ZeroDivisionError:
                    rscu[codon] = 0
                except KeyError:
                    rscu[codon] = 0

        return rscu

    def read_translation_table(self, translation_table: str) -> dict:
        """
        return translation table with codons as keys and amino acids as values
        """

        trans_table = dict()

        if translation_table is None or translation_table == "1":
            pathing = path.join(here, "../../tables/standard_genetic_code.txt")
        elif translation_table == "2":
            pathing = path.join(here, "../../tables/vertebrate_mitochondrial_code.txt")
        elif translation_table == "3":
            pathing = path.join(here, "../../tables/yeast_mitochondrial_code.txt")
        elif translation_table == "4":
            pathing = path.join(
                here,
                "../../tables/mold_protozoan_and_coelenterate_mitochondrial_code_and_the_mycoplasma_spiroplasma.txt",
            )
        elif translation_table == "5":
            pathing = path.join(
                here, "../../tables/invertebrate_mitochondrial_code.txt"
            )
        elif translation_table == "6":
            pathing = path.join(
                here, "../../tables/ciliate_dasycladacean_and_hexamita_nuclear_code.txt"
            )
        elif translation_table == "9":
            pathing = path.join(
                here, "../../tables/echinoderm_and_flatworm_mitochondrial_code.txt"
            )
        elif translation_table == "10":
            pathing = path.join(here, "../../tables/euplotid_nuclear_code.txt")
        elif translation_table == "11":
            pathing = path.join(
                here, "../../tables/bacterial_archaeal_and_plant_plastid_code.txt"
            )
        elif translation_table == "12":
            pathing = path.join(here, "../../tables/alternative_yeast_nuclear_code.txt")
        elif translation_table == "13":
            pathing = path.join(here, "../../tables/ascidian_mitochondrial_code.txt")
        elif translation_table == "14":
            pathing = path.join(
                here, "../../tables/alternative_flatworm_mitochondrial_code.txt"
            )
        elif translation_table == "16":
            pathing = path.join(
                here, "../../tables/chlorophycean_mitochondrial_code.txt"
            )
        elif translation_table == "21":
            pathing = path.join(here, "../../tables/trematode_mitochondrial_code.txt")
        elif translation_table == "22":
            pathing = path.join(
                here, "../../tables/scenedesmus_obliquus_mitochondrial_code.txt"
            )
        elif translation_table == "23":
            pathing = path.join(
                here, "../../tables/thraustochytrium_mitochondrial_code.txt"
            )
        elif translation_table == "24":
            pathing = path.join(
                here, "../../tables/rhabdopleuridae_mitochondrial_code.txt"
            )
        elif translation_table == "25":
            pathing = path.join(
                here, "../../tables/candidate_division_sr1_and_gracilibacteria_code.txt"
            )
        elif translation_table == "26":
            pathing = path.join(
                here, "../../tables/pachysolen_tannophilus_nuclear_code.txt"
            )
        elif translation_table == "27":
            pathing = path.join(here, "../../tables/karyorelict_nuclear_code.txt")
        elif translation_table == "28":
            pathing = path.join(here, "../../tables/condylostoma_nuclear_code.txt")
        elif translation_table == "29":
            pathing = path.join(here, "../../tables/mesodinium_nuclear_code.txt")
        elif translation_table == "30":
            pathing = path.join(here, "../../tables/peritrich_nuclear_code.txt")
        elif translation_table == "31":
            pathing = path.join(here, "../../tables/blastocrithidia_nuclear_code.txt")
        elif translation_table == "33":
            pathing = path.join(
                here, "../../tables/cephalodiscidae_mitochondrial_UAA_tyr_code.txt"
            )
        elif translation_table == "50":
            pathing = path.join(here, "../../tables/CUG_ala_code.txt")
        # case handling for a custom translation table
        else:
            pathing = str(translation_table)

        with open(pathing) as code:
            for line in code:
                line = line.split()
                trans_table[line[0]] = line[1]
        return trans_table
