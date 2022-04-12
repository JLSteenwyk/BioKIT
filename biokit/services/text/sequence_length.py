from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class SequenceLength(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        records = read_and_parse_fasta_seqio(self.fasta)
        for seq_record in records:
            print(f"{seq_record.id}\t{len(seq_record)}")

    def process_args(self, args):
        return dict(fasta=args.fasta)
