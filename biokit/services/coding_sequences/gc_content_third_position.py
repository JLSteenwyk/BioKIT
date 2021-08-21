import re
import sys

from Bio import SeqIO

from .base import CodingSequence


class GCContentThirdPosition(CodingSequence):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = SeqIO.parse(self.fasta, "fasta")

        if self.verbose:
            for record in records:
                third_position_char = []
                third_position_char = self.get_third_position_char(
                    record, third_position_char
                )

                _, matches = self.find_matches("".join(third_position_char))

                print(f"{record.id}\t{round(len(matches)/len(third_position_char), 4)}")

        else:
            third_position_char = []
            number_of_third_positions = 0

            for record in records:
                third_position_char = []
                third_position_char = self.get_third_position_char(
                    record, third_position_char
                )

            _, matches = self.find_matches("".join(third_position_char))

            print(f"{round(len(matches)/len(third_position_char), 4)}")

    def process_args(self, args):
        return dict(fasta=args.fasta, verbose=args.verbose)

    def find_matches(self, seq: str):
        regex_pattern = re.compile("[GgCc]")
        matches = regex_pattern.findall(seq)
        return seq, matches

    def get_third_position_char(self, record, third_position_char: list):
        length_of_coding_seq = len(record._seq)
        for i in range(0, length_of_coding_seq):
            if (i + 1) % 3 == 0:
                third_position_char.append(record._seq[i])
        return third_position_char
