from typing import Any

from Bio import SeqIO

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class RemoveShortSequences(Text):
    def __init__(self, args: Any) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None or self.threshold is None or self.output is None:
            raise ValueError("fasta, threshold, and output cannot be None")
        record_dict = read_and_parse_fasta_seqio(self.fasta)

        out_records: list[Any] = []

        for i in record_dict:
            if len(i.seq) > self.threshold:
                out_records.append(i)

        SeqIO.write(out_records, self.output, "fasta")

    def process_args(self, args: Any) -> dict[str, Any]:
        if args.threshold is None:
            threshold = 500
        else:
            threshold = int(args.threshold)

        if args.output is None:
            output = args.fasta + ".long_seqs.fa"
        else:
            output = args.output

        return dict(
            fasta=args.fasta,
            output=output,
            threshold=threshold
        )
