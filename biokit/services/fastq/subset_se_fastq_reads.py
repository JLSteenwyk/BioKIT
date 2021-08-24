from datetime import datetime
from itertools import islice
import math
import random
import re

from .base import FastQ


class SubsetSEFastQReads(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        reads = []
        with open(self.fastq) as in_handle:
            try:
                while True:
                    read = []
                    data = islice(in_handle, 0, 4)
                    read.append(next(data).strip())
                    read.append(next(data).strip())
                    read.append(next(data).strip())
                    read.append(next(data).strip())

                    reads.append(read)
            except StopIteration:
                1

        number_of_total_reads = len(reads)
        number_of_reads_to_subsample = self.proper_round(
            number_of_total_reads * self.percent
        )
        random.seed(self.seed)
        random_list = random.sample(
            range(0, number_of_total_reads), number_of_reads_to_subsample
        )

        subsetted_reads = []
        for val in random_list:
            for ele in reads[val]:
                subsetted_reads.append(ele)

        # write output file
        with open(self.output_file, "w") as output_fastq_file_name:
            output_fastq_file_name.write("\n".join(subsetted_reads))

    def proper_round(self, num):
        if num - math.floor(num) < 0.5:
            return math.floor(num)
        else:
            return math.ceil(num)

    def process_args(self, args):
        if args.output_file is None:
            output_file = re.sub(".fastq$|.fq$", "_subset.fq", args.fastq)
        else:
            output_file = args.output_file

        if args.percent is None:
            percent = 10 / 100
        else:
            percent = float(args.percent) / 100

        if args.seed is None:
            now = datetime.now()
            now = now.strftime("%H%M%S%m%d%Y")
            seed = int(now)
        else:
            seed = args.seed

        return dict(
            fastq=args.fastq,
            output_file=output_file,
            percent=percent,
            seed=seed,
        )
