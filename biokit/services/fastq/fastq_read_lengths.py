import statistics as stat

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .base import FastQ


class FastQReadLengths(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        read_lens = []

        with open(self.fastq) as in_handle:
            for _, seq, _ in FastqGeneralIterator(in_handle):
                read_lens.append(len(seq))

        if self.verbose:
            for read_len in read_lens:
                print(read_len)
        else:
            mean = round(stat.mean(read_lens), 4)
            stdev = round(stat.stdev(read_lens), 4)
            print(f"{mean} +/- {stdev}")

    def process_args(self, args):
        return dict(
            fastq=args.fastq,
            verbose=args.verbose,
        )
