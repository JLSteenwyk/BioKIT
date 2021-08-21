from Bio import SeqIO

from .base import Text


class ReorderBySequenceLength(Text):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        records = SeqIO.parse(self.fasta, "fasta")

        record_lens = []
        for seq_record in records:
            temp = []
            temp.append(seq_record.id)
            temp.append(len(seq_record))
            record_lens.append(temp)

        record_lens.sort(key=lambda x: x[1], reverse=True)

        records = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta"))

        with open(self.output_file_path, "w") as output_file_path:
            for record in record_lens:
                SeqIO.write(records[record[0]], output_file_path, "fasta")

    def process_args(self, args):
        if args.output is None:
            output_file_path = f"{args.fasta}.reordered.fa"
        else:
            output_file_path = f"{args.output}"

        return dict(fasta=args.fasta, output_file_path=output_file_path)
