import random
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .base import Text
from ...helpers.files import read_and_parse_fasta_seqio


class ShuffleSequences(Text):
    def __init__(self, args: Any) -> None:
        self.seed: int | None = (
            int(args.seed) if getattr(args, "seed", None) is not None else None
        )
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        if self.fasta is None or self.output is None:
            raise ValueError("fasta and output cannot be None")

        if self.seed is not None:
            random.seed(self.seed)

        records = read_and_parse_fasta_seqio(self.fasta)
        shuffled_records: list[SeqRecord] = []

        for record in records:
            chars = list(str(record.seq))
            random.shuffle(chars)
            new_record = SeqRecord(
                Seq("".join(chars)),
                id=record.id,
                description=record.description,
            )
            shuffled_records.append(new_record)

        SeqIO.write(shuffled_records, self.output, "fasta")

    def process_args(self, args: Any) -> dict[str, Any]:
        if args.output is None:
            output = args.fasta + ".shuffled.fa"
        else:
            output = args.output

        return dict(
            fasta=args.fasta,
            output=output,
        )
