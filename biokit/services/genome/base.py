from Bio import SeqIO

from ..base import BaseService


class Genome(BaseService):
    def __init__(
        self,
        *args,
        fasta=None,
        threshold=None,
        verbose=None,
    ):
        self.fasta = fasta
        self.threshold = threshold
        self.verbose = verbose

    def calc_l50(self):
        """
        calculate L50 of a genome assembly
        """
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        contig_lens = []
        for seq_record in records:
            contig_lens.append(len(seq_record))

        # sort and reverse contig lengths
        contig_lens.sort(reverse=True)

        # calculate N50
        sum_contig_lens = sum(contig_lens)
        threshold = sum_contig_lens * 0.50
        n50 = 0
        l50 = 0

        for contig_len in contig_lens:
            n50 += contig_len
            l50 += 1
            if n50 >= threshold:
                return l50
                break

    def calc_l90(self):
        """
        calculate L90 of a genome assembly
        """
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        contig_lens = []
        for seq_record in records:
            contig_lens.append(len(seq_record))

        # sort and reverse contig lengths
        contig_lens.sort(reverse=True)

        # calculate N50
        sum_contig_lens = sum(contig_lens)
        threshold = sum_contig_lens * 0.90
        n90 = 0
        l90 = 0

        for contig_len in contig_lens:
            n90 += contig_len
            l90 += 1
            if n90 >= threshold:
                return l90
                break

    def calc_n50(self):
        """
        calculate n50 of a genome assembly
        """
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        contig_lens = []
        for seq_record in records:
            contig_lens.append(len(seq_record))

        # sort and reverse contig lengths
        contig_lens.sort(reverse=True)

        # calculate N50
        sum_contig_lens = sum(contig_lens)
        n50_threshold = sum_contig_lens * 0.50
        curr = 0

        for contig_len in contig_lens:
            curr += contig_len
            if curr >= n50_threshold:
                return contig_len
                break

    def calc_n90(self):
        """
        calculate n90 of a genome assembly
        """
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        contig_lens = []
        for seq_record in records:
            contig_lens.append(len(seq_record))

        # sort and reverse contig lengths
        contig_lens.sort(reverse=True)

        # calculate N90
        sum_contig_lens = sum(contig_lens)
        threshold = sum_contig_lens * 0.90
        n90 = 0

        for contig_len in contig_lens:
            n90 += contig_len
            if n90 >= threshold:
                return contig_len
                break

    def longest_scaffold(self):
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        max_len = 0
        for seq_record in records:
            if len(seq_record) > max_len:
                max_len = len(seq_record)

        # get longest contig length
        return max_len

    def number_of_large_scaffolds(self):
        """
        calculate the number of large scaffolds
        """
        records = SeqIO.parse(self.fasta, "fasta")
        cnt = 0
        total_len = 0
        for seq_record in records:
            seq_len = len(seq_record)
            if seq_len > self.threshold:
                cnt += 1
                total_len += seq_len
        return cnt, total_len

    def number_of_scaffolds(self):
        """
        count number of scaffolds
        """
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        cnt = 0
        for record in records:
            cnt += 1
        return cnt

    def sum_of_scaffold_lengths(self):
        """
        sum of scaffold lengths
        """
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        cnt = 0
        for record in records:
            cnt += len(record)
        return cnt
