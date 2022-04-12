from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class SequenceComplement(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        records = read_and_parse_fasta_seqio(self.fasta)
        if not self.reverse:
            for seq_record in records:
                print(f">{seq_record.id}\n{seq_record.seq.complement()}")
        else:
            for seq_record in records:
                print(f">{seq_record.id}\n{seq_record.seq.reverse_complement()}")

    def process_args(self, args):
        return dict(
            fasta=args.fasta,
            reverse=args.reverse,
        )
