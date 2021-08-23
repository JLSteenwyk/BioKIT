from datetime import datetime
from itertools import islice
import math
import random
import re

from .base import FastQ


class SubsetPEFastQReads(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        reads1 = []
        reads2 = []
        with open(self.fastq1) as in_handle1, open(self.fastq2) as in_handle2:
            try:
                while True:
                    read1 = []
                    data1 = islice(in_handle1, 0, 4)
                    read1.append(next(data1).strip())
                    read1.append(next(data1).strip())
                    read1.append(next(data1).strip())
                    read1.append(next(data1).strip())

                    reads1.append(read1)

                    read2 = []
                    data2 = islice(in_handle2, 0, 4)
                    read2.append(next(data2).strip())
                    read2.append(next(data2).strip())
                    read2.append(next(data2).strip())
                    read2.append(next(data2).strip())

                    reads2.append(read2)

            except StopIteration:
                1

        number_of_total_reads = len(reads1)
        number_of_reads_to_subsample = self.proper_round(
            number_of_total_reads * self.percent
        )
        random.seed(self.seed)
        random_list = random.sample(
            range(0, number_of_total_reads), number_of_reads_to_subsample
        )

        subsetted_reads1 = []
        subsetted_reads2 = []
        for val in random_list:
            for ele in reads1[val]:
                subsetted_reads1.append(ele)
            for ele in reads2[val]:
                subsetted_reads2.append(ele)

        # write output file
        output_file_1 = re.sub(".fastq$|.fq$", "_subset.fq", self.fastq1)
        with open(output_file_1, "w") as output_fastq_file_name:
            output_fastq_file_name.write("\n".join(subsetted_reads1))
        output_file_2 = re.sub(".fastq$|.fq$", "_subset.fq", self.fastq2)
        with open(output_file_2, "w") as output_fastq_file_name:
            output_fastq_file_name.write("\n".join(subsetted_reads2))

    def proper_round(self, num):
        if num - math.floor(num) < 0.5:
            return math.floor(num)
        else:
            return math.ceil(num)

    def process_args(self, args):
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
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            percent=percent,
            seed=seed,
        )
