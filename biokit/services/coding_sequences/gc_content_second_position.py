import re

from .base import CodingSequence
from ...helpers.files import read_and_parse_fasta_seqio


class GCContentSecondPosition(CodingSequence):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = read_and_parse_fasta_seqio(self.fasta)

        if self.verbose:
            for record in records:
                second_position_char = []
                second_position_char = self.get_second_position_char(
                    record, second_position_char
                )

                _, matches = self.find_matches("".join(second_position_char))

                print(
                    f"{record.id}\t{round(len(matches)/len(second_position_char), 4)}"
                )

        else:
            second_position_char = []

            for record in records:
                second_position_char = []
                second_position_char = self.get_second_position_char(
                    record, second_position_char
                )

            _, matches = self.find_matches("".join(second_position_char))

            print(f"{round(len(matches)/len(second_position_char), 4)}")

    def process_args(self, args):
        return dict(fasta=args.fasta, verbose=args.verbose)

    def find_matches(self, seq: str):
        regex_pattern = re.compile("[GgCc]")
        matches = regex_pattern.findall(seq)
        return seq, matches

    def get_second_position_char(self, record, second_position_char: list):
        length_of_coding_seq = len(record._seq)
        for i in range(0, length_of_coding_seq, 3):
            try:
                second_position_char.append(record._seq[i:i + 3][1])
            except IndexError:
                continue
        return second_position_char
