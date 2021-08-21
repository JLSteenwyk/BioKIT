from os import path

from ..base import BaseService

here = path.dirname(__file__)


class Text(BaseService):
    def __init__(
        self,
        *args,
        ambiguous_character=None,
        entry=None,
        fasta=None,
        idmap=None,
        input_file=None,
        input_file_format=None,
        output=None,
        output_file=None,
        output_file_format=None,
        output_file_path=None,
        reverse=None,
        threshold=None,
    ):
        self.ambiguous_character = ambiguous_character
        self.entry = entry
        self.fasta = fasta
        self.idmap = idmap
        self.input_file = input_file
        self.input_file_format = input_file_format
        self.output = output
        self.output_file = output_file
        self.output_file_format = output_file_format
        self.output_file_path = output_file_path
        self.reverse = reverse
        self.threshold = threshold

    def read_translation_table(self) -> dict:
        """
        return translation table with codons as keys and amino acids as values
        """

        trans_table = dict()

        if self.translation_table is None or self.translation_table == "1":
            pathing = path.join(here, "../../tables/standard_genetic_code.txt")
        elif self.translation_table == "2":
            pathing = path.join(here, "../../tables/vertebrate_mitochondrial_code.txt")
        elif self.translation_table == "3":
            pathing = path.join(here, "../../tables/yeast_mitochondrial_code.txt")
        elif self.translation_table == "4":
            pathing = path.join(
                here,
                "../../tables/mold_protozoan_and_coelenterate_mitochondrial_code_and_the_mycoplasma_spiroplasma.txt",
            )
        elif self.translation_table == "5":
            pathing = path.join(
                here, "../../tables/invertebrate_mitochondrial_code.txt"
            )
        elif self.translation_table == "6":
            pathing = path.join(
                here, "../../tables/ciliate_dasycladacean_and_hexamita_nuclear_code.txt"
            )
        elif self.translation_table == "9":
            pathing = path.join(
                here, "../../tables/echinoderm_and_flatworm_mitochondrial_code.txt"
            )
        elif self.translation_table == "10":
            pathing = path.join(here, "../../tables/euplotid_nuclear_code.txt")
        elif self.translation_table == "11":
            pathing = path.join(
                here, "../../tables/bacterial_archaeal_and_plant_plastid_code.txt"
            )
        elif self.translation_table == "12":
            pathing = path.join(here, "../../tables/alternative_yeast_nuclear_code.txt")
        elif self.translation_table == "13":
            pathing = path.join(here, "../../tables/ascidian_mitochondrial_code.txt")
        elif self.translation_table == "14":
            pathing = path.join(
                here, "../../tables/alternative_flatworm_mitochondrial_code.txt"
            )
        elif self.translation_table == "16":
            pathing = path.join(
                here, "../../tables/chlorophycean_mitochondrial_code.txt"
            )
        elif self.translation_table == "21":
            pathing = path.join(here, "../../tables/trematode_mitochondrial_code.txt")
        elif self.translation_table == "22":
            pathing = path.join(
                here, "../../tables/scenedesmus_obliquus_mitochondrial_code.txt"
            )
        elif self.translation_table == "23":
            pathing = path.join(
                here, "../../tables/thraustochytrium_mitochondrial_code.txt"
            )
        elif self.translation_table == "24":
            pathing = path.join(
                here, "../../tables/rhabdopleuridae_mitochondrial_code.txt"
            )
        elif self.translation_table == "25":
            pathing = path.join(
                here, "../../tables/candidate_division_sr1_and_gracilibacteria_code.txt"
            )
        elif self.translation_table == "26":
            pathing = path.join(
                here, "../../tables/pachysolen_tannophilus_nuclear_code.txt"
            )
        elif self.translation_table == "27":
            pathing = path.join(here, "../../tables/karyorelict_nuclear_code.txt")
        elif self.translation_table == "28":
            pathing = path.join(here, "../../tables/condylostoma_nuclear_code.txt")
        elif self.translation_table == "29":
            pathing = path.join(here, "../../tables/mesodinium_nuclear_code.txt")
        elif self.translation_table == "30":
            pathing = path.join(here, "../../tables/peritrich_nuclear_code.txt")
        elif self.translation_table == "31":
            pathing = path.join(here, "../../tables/blastocrithidia_nuclear_code.txt")
        elif self.translation_table == "33":
            pathing = path.join(
                here, "../../tables/cephalodiscidae_mitochondrial_UAA_tyr_code.txt"
            )
        elif self.translation_table == "50":
            pathing = path.join(here, "../../tables/CUG_ala_code.txt")
        # case handling for a custom translation table
        else:
            pathing = str(self.translation_table)

        with open(pathing) as code:
            for line in code:
                line = line.split()
                trans_table[line[0]] = line[1]
        return trans_table
