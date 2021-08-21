import statistics as stat

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .base import FastQ


class FastQReadLengths(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        read_lens = []
        with open(self.fastq) as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                read_lens.append(len(seq))

        if self.verbose:
            try:
                for read_len in read_lens:
                    print(read_len)
            except BrokenPipeError:
                pass
        else:
            mean = stat.mean(read_lens)
            stdev = stat.stdev(read_lens)
            print(f"{mean} +/- {stdev}")

    def process_args(self, args):
        return dict(
            fastq=args.fastq,
            verbose=args.verbose,
        )
