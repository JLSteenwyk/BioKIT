from ..base import BaseService
from ...helpers.files import iter_fasta_sequences


class Genome(BaseService):
    def __init__(
        self,
        *args,
        fasta=None,
        threshold=None,
        verbose=None,
        **kwargs,
    ):
        self.fasta = fasta
        self.threshold = threshold
        self.verbose = verbose
        self.output_format = kwargs.get("output_format")
        self._contig_lengths = None
        self._base_counts = None

    def _ensure_cached_genome_stats(self):
        if self._contig_lengths is None or self._base_counts is None:
            contig_lengths = []
            base_counts = {"A": 0, "T": 0, "G": 0, "C": 0}
            for seq in iter_fasta_sequences(self.fasta):
                seq = seq.upper()
                contig_lengths.append(len(seq))
                base_counts["A"] += seq.count("A")
                base_counts["T"] += seq.count("T")
                base_counts["G"] += seq.count("G")
                base_counts["C"] += seq.count("C")
            self._contig_lengths = contig_lengths
            self._base_counts = base_counts

    def _get_contig_lengths(self):
        self._ensure_cached_genome_stats()
        return self._contig_lengths

    def _get_base_counts(self):
        self._ensure_cached_genome_stats()
        return self._base_counts

    def calc_l50(self):
        """
        calculate L50 of a genome assembly
        """
        contig_lens = self._get_contig_lengths()[:]

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

    def calc_l90(self):
        """
        calculate L90 of a genome assembly
        """
        contig_lens = self._get_contig_lengths()[:]

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

    def calc_n50(self):
        """
        calculate n50 of a genome assembly
        """
        contig_lens = self._get_contig_lengths()[:]

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

    def calc_n90(self):
        """
        calculate n90 of a genome assembly
        """
        contig_lens = self._get_contig_lengths()[:]

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

    def longest_scaffold(self):
        contig_lens = self._get_contig_lengths()
        return max(contig_lens, default=0)

    def number_of_large_scaffolds(self):
        """
        calculate the number of large scaffolds
        """
        cnt = 0
        total_len = 0
        for seq_len in self._get_contig_lengths():
            if seq_len > self.threshold:
                cnt += 1
                total_len += seq_len
        return cnt, total_len

    def number_of_scaffolds(self):
        """
        count number of scaffolds
        """
        return len(self._get_contig_lengths())

    def sum_of_scaffold_lengths(self):
        """
        sum of scaffold lengths
        """
        return sum(self._get_contig_lengths())
