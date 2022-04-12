from os import path
import re

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .base import FastQ

here = path.dirname(__file__)


class TrimSEAdaptersFastQ(FastQ):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        adapter_table = self.read_adapter_table(self.adapters)  # noqa

        good_reads = []
        cnt = 0
        kept = 0
        removed = 0
        adapter_removed = 0

        with open(self.fastq) as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                seq_record_len = len(seq)  # cache this for later
                if seq_record_len < self.length:
                    # Too short to keep
                    continue
                for _, adapter_seq in adapter_table.items():
                    try:
                        if index < seq.find(adapter_seq) and index != -1: # noqa
                            index = seq.find(adapter_seq)
                            len_adaptor = len(adapter_seq)
                    except UnboundLocalError:
                        index = seq.find(adapter_seq)
                        len_adaptor = len(adapter_seq)

                if index == -1:
                    # adaptor not found, so won't trim
                    good_reads.append("@" + title)
                    good_reads.append(seq)
                    good_reads.append("+" + title)
                    good_reads.append(qual)
                    kept += 1
                    cnt += 1
                if index == 5000:
                    index = 0
                elif seq_record_len - index - len_adaptor >= self.length:
                    # after trimming this will still be long enough
                    good_reads.append("@" + title)
                    good_reads.append(seq[index + len_adaptor:])
                    good_reads.append("+" + title)
                    good_reads.append(qual)
                    kept += 1
                    cnt += 1
                    adapter_removed += 1

        print(
            f"Reads processed: {cnt}\nReads kept: {kept}\nReads removed: {removed}\nAdapaters removed: {adapter_removed}"
        )

        # write output file
        with open(self.output_file, "w") as output_fastq_file_name:
            output_fastq_file_name.write("\n".join(good_reads))

    def process_args(self, args):
        if args.adapters is None:
            adapters = "TruSeq2-SE"
        else:
            adapters = args.adapters

        if args.length is None:
            length = 20
        else:
            length = int(args.length)

        if args.output_file is None:
            output_file = re.sub(".fastq$|.fq$", "_adapter_trimmed.fq", args.fastq)
        else:
            output_file = args.output_file

        return dict(
            fastq=args.fastq,
            adapters=adapters,
            length=length,
            output_file=output_file,
        )
