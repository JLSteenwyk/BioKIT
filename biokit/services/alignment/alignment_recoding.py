from os import path
import sys

from .base import Alignment
from ...helpers.files import read_alignment_alignio

here = path.dirname(__file__)


class AlignmentRecoding(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment = read_alignment_alignio(self.fasta)

        recoding_table = self.read_recoding_table(self.code)

        recoded_alignment = self.recode_alignment_as_dict(alignment, recoding_table)

        for k, v in recoded_alignment.items():
            print(f">{k}\n{''.join(v)}")

    def recode_alignment_as_dict(self, alignment, recoding_table: dict) -> dict:
        recoded_alignment = dict()
        for i in range(0, len(alignment)):
            recoded_sequence_i = []
            for j in range(alignment.get_alignment_length()):
                sequence_ij = alignment[i, j].upper()
                if sequence_ij in ["?", "-"]:
                    recoded_sequence_i.append(sequence_ij)
                else:
                    recoded_sequence_i.append(recoding_table[sequence_ij])

            recoded_alignment[alignment[i].id] = recoded_sequence_i

        return(recoded_alignment)

    def read_recoding_table(self, recoding: str) -> dict:
        """
        return translation table with codons as keys and amino acids as values
        """

        recoding_table = dict()

        if recoding is None:
            print("Please specify a recoding table")
            sys.exit()
        elif recoding == "RY-nucleotide":
            pathing = path.join(here, "../../recoding_tables/RY-nucleotide.txt")
        elif recoding == "Dayhoff-6":
            pathing = path.join(here, "../../recoding_tables/Dayhoff-6.txt")
        elif recoding == "SandR-6":
            pathing = path.join(here, "../../recoding_tables/S_and_R-6.txt")
        elif recoding == "KGB-6":
            pathing = path.join(here, "../../recoding_tables/KGB-6.txt")
        # case handling for a custom translation table
        else:
            pathing = str(recoding)

        with open(pathing) as code:
            for line in code:
                line = line.split()
                if line[1].upper() in recoding_table.keys():
                    recoding_table[line[1]].upper().append(line[0].upper())
                else:
                    recoding_table[line[1]] = line[0].upper()

        return recoding_table

    def process_args(self, args):
        return dict(
            fasta=args.fasta,
            code=args.code
        )
