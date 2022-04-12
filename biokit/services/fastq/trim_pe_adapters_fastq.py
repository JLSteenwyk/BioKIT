from os import path
import re

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .base import FastQ

here = path.dirname(__file__)


class TrimPEAdaptersFastQ(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        adapter_table = self.read_adapter_table(self.adapters)  # noqa

        good_reads_paired_1 = []
        good_reads_paired_2 = []
        good_reads_unpaired_1 = []
        good_reads_unpaired_2 = []
        cnt = 0
        kept = 0
        removed = 0
        adapter_removed = 0

        with open(self.fastq1) as in_handle1, open(self.fastq2) as in_handle2:
            for (title1, seq1, qual1), (title2, seq2, qual2) in zip(
                FastqGeneralIterator(in_handle1), FastqGeneralIterator(in_handle2)
            ):

                # logic for keeping a read or not
                keep_1 = True
                keep_2 = True

                seq_record_len = len(seq1)  # cache this for later
                if seq_record_len < self.length:
                    # Too short to keep
                    keep_1 = False
                for _, adapter_seq in adapter_table.items():
                    try:
                        if index1 < seq1.find(adapter_seq) and index1 != -1: # noqa
                            index1 = seq1.find(adapter_seq)
                            len_adaptor = len(adapter_seq)
                    except UnboundLocalError:
                        index1 = seq1.find(adapter_seq)
                        len_adaptor = len(adapter_seq)

                if seq_record_len - index1 - len_adaptor < self.length:
                    keep_1 = False

                seq_record_len = len(seq2)  # cache this for later
                if seq_record_len < self.length:
                    # Too short to keep
                    keep_2 = False
                for _, adapter_seq in adapter_table.items():
                    try:
                        if index2 < seq2.find(adapter_seq) and index1 != -1: # noqa
                            index2 = seq2.find(adapter_seq)
                            len_adaptor = len(adapter_seq)
                    except UnboundLocalError:
                        index2 = seq2.find(adapter_seq)
                        len_adaptor = len(adapter_seq)

                if seq_record_len - index2 - len_adaptor < self.length:
                    keep_2 = False

                # keep both reads
                if keep_1 and keep_2:
                    if index1 == -1:
                        # adaptor not found, so won't trim
                        good_reads_paired_1.append("@" + title1)
                        good_reads_paired_1.append(seq1)
                        good_reads_paired_1.append("+" + title1)
                        good_reads_paired_1.append(qual1)
                        kept += 1

                    elif seq_record_len - index1 - len_adaptor >= self.length:
                        # after trimming this will still be long enough
                        good_reads_paired_1.append("@" + title1)
                        good_reads_paired_1.append(seq1[index1 + len_adaptor:])
                        good_reads_paired_1.append("+" + title1)
                        good_reads_paired_1.append(qual1)
                        kept += 1
                        adapter_removed += 1

                    if index2 == -1:
                        # adaptor not found, so won't trim
                        good_reads_paired_2.append("@" + title2)
                        good_reads_paired_2.append(seq2)
                        good_reads_paired_2.append("+" + title2)
                        good_reads_paired_2.append(qual2)
                        kept += 1

                    elif seq_record_len - index2 - len_adaptor >= self.length:
                        # after trimming this will still be long enough
                        good_reads_paired_2.append("@" + title2)
                        good_reads_paired_2.append(seq2[index2 + len_adaptor:])
                        good_reads_paired_2.append("+" + title2)
                        good_reads_paired_2.append(qual2)
                        kept += 1
                        adapter_removed += 1

                # keep 1 but not 2
                elif keep_1 and not keep_2:
                    if index1 == -1:
                        # adaptor not found, so won't trim
                        good_reads_unpaired_1.append("@" + title1)
                        good_reads_unpaired_1.append(seq1)
                        good_reads_unpaired_1.append("+" + title1)
                        good_reads_unpaired_1.append(qual1)
                        kept += 1

                    elif seq_record_len - index1 - len_adaptor >= self.length:
                        # after trimming this will still be long enough
                        good_reads_unpaired_1.append("@" + title1)
                        good_reads_unpaired_1.append(seq1[index1 + len_adaptor:])
                        good_reads_unpaired_1.append("+" + title1)
                        good_reads_unpaired_1.append(qual1)
                        kept += 1
                        adapter_removed += 1
                # keep 2 but not 1
                else:
                    if index2 == -1:
                        # adaptor not found, so won't trim
                        good_reads_unpaired_2.append("@" + title2)
                        good_reads_unpaired_2.append(seq2)
                        good_reads_unpaired_2.append("+" + title2)
                        good_reads_unpaired_2.append(qual2)
                        kept += 1

                    elif seq_record_len - index2 - len_adaptor >= self.length:
                        # after trimming this will still be long enough
                        good_reads_unpaired_2.append("@" + title2)
                        good_reads_unpaired_2.append(seq2[index2 + len_adaptor:])
                        good_reads_unpaired_2.append("+" + title2)
                        good_reads_unpaired_2.append(qual2)
                        kept += 1
                        adapter_removed += 1
                cnt += 1
        print(
            f"Reads processed: {cnt}\nReads kept: {kept}\nReads removed: {removed}\nAdapaters removed: {adapter_removed}"
        )

        # write output files
        # write paired output files
        if len(good_reads_paired_1) != 0:
            output_file_paired_1 = re.sub(
                ".fastq$|.fq$", "_paired_adapter_trimmed.fq", self.fastq1
            )
            with open(output_file_paired_1, "w") as output_fastq_file_name_1:
                output_fastq_file_name_1.write("\n".join(good_reads_paired_1))
        if len(good_reads_paired_2) != 0:
            output_file_paired_2 = re.sub(
                ".fastq$|.fq$", "_paired_adapter_trimmed.fq", self.fastq2
            )
            with open(output_file_paired_2, "w") as output_fastq_file_name_2:
                output_fastq_file_name_2.write("\n".join(good_reads_paired_2))

        # write unpaired output files
        if len(good_reads_unpaired_1) != 0:
            output_file_unpaired_1 = re.sub(
                ".fastq$|.fq$", "_unpaired_adapter_trimmed.fq", self.fastq1
            )
            with open(output_file_unpaired_1, "w") as output_fastq_file_name_1:
                output_fastq_file_name_1.write("\n".join(good_reads_unpaired_1))
        if len(good_reads_unpaired_2) != 0:
            output_file_unpaired_2 = re.sub(
                ".fastq$|.fq$", "_unpaired_adapter_trimmed.fq", self.fastq2
            )
            with open(output_file_unpaired_2, "w") as output_fastq_file_name_2:
                output_fastq_file_name_2.write("\n".join(good_reads_unpaired_2))

    def process_args(self, args):
        if args.adapters is None:
            adapters = "TruSeq2-PE"
        else:
            adapters = args.adapters

        if args.length is None:
            length = 20
        else:
            length = int(args.length)

        return dict(
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            adapters=adapters,
            length=length,
        )
