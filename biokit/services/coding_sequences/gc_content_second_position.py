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

                gc_content = self.calculate_gc_content("".join(second_position_char))
                print(f"{record.id}\t{gc_content}")

        else:
            second_position_char = []

            for record in records:
                second_position_char = []
                second_position_char = self.get_second_position_char(
                    record, second_position_char
                )

            gc_content = self.calculate_gc_content("".join(second_position_char))
            print(f"{gc_content}")

    def process_args(self, args):
        return dict(fasta=args.fasta, verbose=args.verbose)

    def calculate_gc_content(self, seq: str) -> float:
        if not seq:
            return 0
        gc_count = seq.count("G") + seq.count("g") + seq.count("C") + seq.count("c")
        return round(gc_count / len(seq), 4)

    def get_second_position_char(self, record, second_position_char: list):
        sequence = str(record.seq)
        second_position_char.extend(sequence[1::3])
        return second_position_char
