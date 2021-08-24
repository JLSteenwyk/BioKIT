from os import path
import re

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .base import FastQ

here = path.dirname(__file__)


class TrimPEFastQ(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # read in quality table
        quality_table = dict()
        pathing = path.join(here, "../../tables/ascii_base_33.txt")
        with open(pathing) as asci_base_33:
            for line in asci_base_33:
                line = line.split()
                quality_table[line[0]] = line[1]

        # initialize lists for keeping outputted reads
        good_reads_paired1 = []
        good_reads_paired2 = []
        good_reads_unpaired1 = []
        good_reads_unpaired2 = []

        # variables for keeping track of counts among processed reads
        cnt = 0
        kept = 0
        removed = 0
        pairs_kept = 0

        # loop through both fastq files
        with open(self.fastq1) as in_handle1, open(self.fastq2) as in_handle2:
            for (title1, seq1, qual1), (title2, seq2, qual2) in zip(
                FastqGeneralIterator(in_handle1), FastqGeneralIterator(in_handle2)
            ):
                # logic for keeping a read or not
                keep_1 = True
                keep_2 = True

                # indexes for where you trim
                trim_idx1 = None
                trim_idx2 = None

                ## read 1 handling # noqa
                # determine index to trim read at
                for i in range(len(seq1)):
                    base_quality = quality_table[qual1[i]]
                    if int(base_quality) < self.minimum:
                        trim_idx1 = i
                        break
                # if there is no base below the quality threshold
                # save the whole read
                if not trim_idx1 and trim_idx1 != 0:
                    kept += 1
                    cnt += 1
                # if the trimming idx is longer than the length
                # threshold, trim the read
                elif trim_idx1 >= self.length:
                    kept += 1
                    cnt += 1
                else:
                    removed += 1
                    cnt += 1
                    keep_1 = False

                ## read 2 handling # noqa
                for j in range(len(seq2)):
                    base_quality = quality_table[qual2[j]]
                    if int(base_quality) < self.minimum:
                        trim_idx2 = j
                        break
                # if there is no base below the quality threshold
                # save the whole read
                if not trim_idx2 and trim_idx2 != 0:
                    kept += 1
                    cnt += 1
                # if the trimming idx is longer than the length
                # threshold, trim the read
                elif trim_idx2 >= self.length:
                    kept += 1
                    cnt += 1
                else:
                    removed += 1
                    cnt += 1
                    keep_2 = False

                # save good reads to lists
                if keep_1 and keep_2:
                    seq1 = seq1[:i]
                    qual1 = qual1[:i]
                    good_reads_paired1.append("@" + title1)
                    good_reads_paired1.append(seq1)
                    good_reads_paired1.append("+" + title1)
                    good_reads_paired1.append(qual1)

                    seq2 = seq2[:j]
                    qual2 = qual2[:j]
                    good_reads_paired2.append("@" + title2)
                    good_reads_paired2.append(seq2)
                    good_reads_paired2.append("+" + title2)
                    good_reads_paired2.append(qual2)

                    kept += 1
                    pairs_kept += 1

                elif keep_1 and not keep_2:
                    seq1 = seq1[:i]
                    qual1 = qual1[:i]
                    good_reads_unpaired1.append("@" + title1)
                    good_reads_unpaired1.append(seq1)
                    good_reads_unpaired1.append("+" + title1)
                    good_reads_unpaired1.append(qual1)

                # logic == not keep_1 and keep_2
                else:
                    seq2 = seq2[:j]
                    qual2 = qual2[:j]
                    good_reads_unpaired2.append("@" + title2)
                    good_reads_unpaired2.append(seq2)
                    good_reads_unpaired2.append("+" + title2)
                    good_reads_unpaired2.append(qual2)

        print(
            f"Reads processed: {cnt}\nReads kept: {kept}\nReads removed: {removed}\n\nPairs kept: {pairs_kept}"
        )

        # write output files
        # write paired output files
        if len(good_reads_paired1) != 0:
            output_file_paired_1 = re.sub(
                ".fastq$|.fq$", "_paired_trimmed.fq", self.fastq1
            )
            with open(output_file_paired_1, "w") as output_fastq_file_name_1:
                output_fastq_file_name_1.write("\n".join(good_reads_paired1))
        if len(good_reads_paired2) != 0:
            output_file_paired_2 = re.sub(
                ".fastq$|.fq$", "_paired_trimmed.fq", self.fastq2
            )
            with open(output_file_paired_2, "w") as output_fastq_file_name_2:
                output_fastq_file_name_2.write("\n".join(good_reads_paired2))

        # write unpaired output files
        if len(good_reads_unpaired1) != 0:
            output_file_unpaired_1 = re.sub(
                ".fastq$|.fq$", "_unpaired_trimmed.fq", self.fastq1
            )
            with open(output_file_unpaired_1, "w") as output_fastq_file_name_1:
                output_fastq_file_name_1.write("\n".join(good_reads_unpaired1))
        if len(good_reads_unpaired2) != 0:
            output_file_unpaired_2 = re.sub(
                ".fastq$|.fq$", "_unpaired_trimmed.fq", self.fastq2
            )
            with open(output_file_unpaired_2, "w") as output_fastq_file_name_2:
                output_fastq_file_name_2.write("\n".join(good_reads_unpaired2))

    def process_args(self, args):
        if args.minimum is None:
            minimum = 20
        else:
            minimum = int(args.minimum)

        if args.length is None:
            length = 20
        else:
            length = int(args.length)

        return dict(
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            minimum=minimum,
            length=length,
        )
