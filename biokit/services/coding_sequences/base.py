from os import path
from collections import Counter
from typing import Any

from ..base import BaseService
from ...helpers.files import iter_fasta_entries

here = path.dirname(__file__)


class CodingSequence(BaseService):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.fasta: str | None = kwargs.get("fasta")
        self.verbose: bool | None = kwargs.get("verbose")
        self.translation_table: str | None = kwargs.get("translation_table")
        self.output_file_path: str | None = kwargs.get("output_file_path")
        self.output_format: str | None = kwargs.get("output_format")

    def calculate_rscu(self, translation_table: dict[str, str]) -> dict[str, float | int]:
        # get codon_table
        codon_table: dict[str, list[str]] = {}
        for k, v in translation_table.items():
            if v in codon_table.keys():
                codon_table[v].append(k)
            else:
                codon_table[v] = [k]

        # get counts of codons
        codon_counts: Counter[str] = Counter()
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        valid_codons = set(translation_table.keys())
        for _, seq in iter_fasta_entries(self.fasta):
            sequence = seq.upper().replace("T", "U")
            if len(sequence) % 3 != 0:
                continue
            for position in range(0, len(sequence), 3):
                codon = sequence[position:position + 3]
                if codon in valid_codons:
                    codon_counts[codon] += 1

        # calculate rscu
        rscu: dict[str, float | int] = {}
        for _, codons in codon_table.items():
            observed_sum = 0
            for codon in codons:
                observed_sum += codon_counts.get(codon, 0)
            for codon in codons:
                codon_count = codon_counts.get(codon, 0)
                if codon_count == 0:
                    rscu[codon] = 0
                    continue
                try:
                    rscu[codon] = round(codon_count / (observed_sum / len(codons)), 4)
                except ZeroDivisionError:
                    rscu[codon] = 0

        return rscu

    @staticmethod
    def calculate_gc_content(seq: str) -> float | int:
        if not seq:
            return 0
        gc_count = seq.count("G") + seq.count("g") + seq.count("C") + seq.count("c")
        return round(gc_count / len(seq), 4)

    @staticmethod
    def get_codon_position_chars(sequence: str, codon_position: int) -> str:
        return sequence[codon_position::3]

    def gc_content_by_codon_position(
        self, codon_position: int
    ) -> list[str]:
        if self.fasta is None:
            raise ValueError("fasta cannot be None")
        output_lines: list[str] = []
        aggregated_chars: list[str] = []

        for seq_id, seq in iter_fasta_entries(self.fasta):
            position_chars = self.get_codon_position_chars(seq, codon_position)
            if self.verbose:
                gc_content = self.calculate_gc_content(position_chars)
                output_lines.append(f"{seq_id}\t{gc_content}")
            else:
                aggregated_chars.extend(position_chars)

        if not self.verbose:
            output_lines.append(str(self.calculate_gc_content("".join(aggregated_chars))))

        return output_lines

    def read_translation_table(self, translation_table: str | None) -> dict[str, str]:
        """
        return translation table with codons as keys and amino acids as values
        """

        trans_table: dict[str, str] = {}

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
                parts = line.split()
                trans_table[parts[0]] = parts[1]
        return trans_table
