from ..base import BaseService


class FastQ(BaseService):
    def __init__(
        self,
        *args,
        fastq=None,
        fastq1=None,
        fastq2=None,
        length=None,
        minimum=None,
        output_file=None,
        percent=None,
        seed=None,
        verbose=None,
    ):
        self.fastq = fastq
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.length = length
        self.minimum = minimum
        self.output_file = output_file
        self.percent = percent
        self.seed = seed
        self.verbose = verbose
