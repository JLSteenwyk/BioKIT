import statistics as stat
import sys

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .base import FastQ


class FastQReadLengths(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        output_format = self.normalize_output_format(self.output_format)
        read_lens = []

        if self.fastq == '-':
            for _, seq, _ in FastqGeneralIterator(sys.stdin):
                read_lens.append(len(seq))
        else:
            with open(self.fastq) as in_handle:
                for _, seq, _ in FastqGeneralIterator(in_handle):
                    read_lens.append(len(seq))

        if self.verbose:
            if output_format == "tsv":
                for read_len in read_lens:
                    print(read_len)
            else:
                rows = [
                    {"read_index": index + 1, "length": read_len}
                    for index, read_len in enumerate(read_lens)
                ]
                print(self.format_rows(rows, output_format))
        else:
            mean = round(stat.mean(read_lens), 4)
            stdev = round(stat.stdev(read_lens), 4)
            if output_format == "tsv":
                print(f"{mean} +/- {stdev}")
            else:
                print(self.format_object({"mean": mean, "stdev": stdev}, output_format))

    def process_args(self, args):
        return dict(
            fastq=args.fastq,
            verbose=args.verbose,
            output_format=getattr(args, "format", None),
        )
