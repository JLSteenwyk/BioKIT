from .base import CodingSequence
from ...helpers.files import read_and_parse_fasta_seqio


class GCContentFirstPosition(CodingSequence):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = read_and_parse_fasta_seqio(self.fasta)

        if self.verbose:
            for record in records:
                first_position_char = []
                first_position_char = self.get_first_position_char(
                    record, first_position_char
                )

                gc_content = self.calculate_gc_content("".join(first_position_char))
                print(f"{record.id}\t{gc_content}")

        else:
            first_position_char = []

            for record in records:
                first_position_char = []
                first_position_char = self.get_first_position_char(
                    record, first_position_char
                )

            gc_content = self.calculate_gc_content("".join(first_position_char))
            print(f"{gc_content}")

    def process_args(self, args):
        return dict(fasta=args.fasta, verbose=args.verbose)

    def calculate_gc_content(self, seq: str) -> float:
        if not seq:
            return 0
        gc_count = seq.count("G") + seq.count("g") + seq.count("C") + seq.count("c")
        return round(gc_count / len(seq), 4)

    def get_first_position_char(self, record, first_position_char: list):
        sequence = str(record.seq)
        first_position_char.extend(sequence[::3])
        return first_position_char
