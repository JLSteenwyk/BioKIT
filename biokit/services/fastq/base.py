from os import path

from ..base import BaseService

here = path.dirname(__file__)


class FastQ(BaseService):
    def __init__(
        self,
        *args,
        adapters=None,
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
        self.adapters = adapters
        self.fastq = fastq
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.length = length
        self.minimum = minimum
        self.output_file = output_file
        self.percent = percent
        self.seed = seed
        self.verbose = verbose

    def read_adapter_table(self, adapters: str) -> dict:
        """
        return adapter table with adapter name as key and adapter sequence as value
        """

        adapter_table = dict()

        if adapters is None or adapters == "TruSeq2-SE":
            pathing = path.join(here, "../../adapters/TruSeq2-SE.txt")
        elif adapters == "TruSeq2-PE":
            pathing = path.join(here, "../../adapters/TruSeq2-PE.txt")
        elif adapters == "TruSeq3-PE-2":
            pathing = path.join(here, "../../adapters/TruSeq3-PE-2.txt")
        elif adapters == "TruSeq3-PE":
            pathing = path.join(here, "../../adapters/TruSeq3-PE.txt")
        elif adapters == "TruSeq3-SE":
            pathing = path.join(here, "../../adapters/TruSeq3-SE.txt")
        elif adapters == "NexteraPE-PE":
            pathing = path.join(here, "../../adapters/NexteraPE-PE.txt")
        # case handling for a custom translation table
        else:
            pathing = str(adapters)

        with open(pathing) as code:
            for line in code:
                line = line.split()
                adapter_table[line[0]] = line[1]
        return adapter_table
