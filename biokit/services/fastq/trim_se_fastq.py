from os import path
import re

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .base import FastQ

here = path.dirname(__file__)


class TrimSEFastQ(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        quality_table = dict()
        pathing = path.join(here, "../../tables/ascii_base_33.txt")
        with open(pathing) as asci_base_33:
            for line in asci_base_33:
                line = line.split()
                quality_table[line[0]] = line[1]

        good_reads = []
        cnt = 0
        kept = 0
        removed = 0
        with open(self.fastq) as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                trim_idx = None
                # determine index to trim read at
                for i in range(len(seq)):
                    base_quality = quality_table[qual[i]]
                    if int(base_quality) < self.minimum:
                        trim_idx = i
                        break
                # if there is no base below the quality threshold
                # save the whole read
                if not trim_idx and trim_idx != 0:
                    good_reads.append("@" + title)
                    good_reads.append(seq)
                    good_reads.append("+" + title)
                    good_reads.append(qual)
                    kept += 1
                    cnt += 1
                # if the trimming idx is longer than the length
                # threshold, trim the read
                elif trim_idx >= self.length:
                    seq = seq[:i]
                    qual = qual[:i]
                    good_reads.append("@" + title)
                    good_reads.append(seq)
                    good_reads.append("+" + title)
                    good_reads.append(qual)
                    kept += 1
                    cnt += 1
                else:
                    removed += 1
                    cnt += 1

        print(f"Reads processed: {cnt}\nReads kept: {kept}\nReads removed: {removed}")

        # write output file
        with open(self.output_file, "w") as output_fastq_file_name:
            output_fastq_file_name.write("\n".join(good_reads))

    def process_args(self, args):
        if args.minimum is None:
            minimum = 20
        else:
            minimum = int(args.minimum)

        if args.length is None:
            length = 20
        else:
            length = int(args.length)

        if args.output_file is None:
            output_file = re.sub(".fastq$|.fq$", "_trimmed.fq", args.fastq)
        else:
            output_file = args.output_file

        return dict(
            fastq=args.fastq,
            minimum=minimum,
            length=length,
            output_file=output_file,
        )
