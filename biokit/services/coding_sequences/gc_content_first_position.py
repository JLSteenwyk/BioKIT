import re

from Bio import SeqIO

from .base import CodingSequence


class GCContentFirstPosition(CodingSequence):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = SeqIO.parse(self.fasta, "fasta")

        if self.verbose:
            for record in records:
                first_position_char = []
                first_position_char = self.get_first_position_char(
                    record, first_position_char
                )

                _, matches = self.find_matches("".join(first_position_char))

                print(f"{record.id}\t{round(len(matches)/len(first_position_char), 4)}")

        else:
            first_position_char = []

            for record in records:
                first_position_char = []
                first_position_char = self.get_first_position_char(
                    record, first_position_char
                )

            _, matches = self.find_matches("".join(first_position_char))

            print(f"{round(len(matches)/len(first_position_char), 4)}")

    def process_args(self, args):
        return dict(fasta=args.fasta, verbose=args.verbose)

    def find_matches(self, seq: str):
        regex_pattern = re.compile("[GgCc]")
        matches = regex_pattern.findall(seq)
        return seq, matches

    def get_first_position_char(self, record, first_position_char: list):
        length_of_coding_seq = len(record._seq)
        for i in range(0, length_of_coding_seq, 3):
            first_position_char.append(record._seq[i : i + 3][0])
        return first_position_char
