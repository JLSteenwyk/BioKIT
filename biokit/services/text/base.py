from collections import Counter
import sys

from Bio import SeqIO

from ..base import BaseService
from ...helpers.files import get_alignment_and_format as get_alignment_and_format_helper

class Text(BaseService):
    def __init__(
        self,
        *args,
        alignment_file_path=None,
        fasta=None,
        output_file_path=None,
        protein_file_path=None,
        nucleotide_file_path=None,
        alignment_list_path=None,
        prefix=None,
        idmap=None,
        reference=None,
        verbose=None,
        entry=None,
        input_file=None,
        output_file_format=None,
        input_file_format=None,
        output_file=None,
        threshold=None,
        ambiguous_character=None,
        reverse=None,
        translation_table=None,
        fastq=None,
        minimum=None,
        length=None,
        output_file_1=None,
        output_file_2=None,
        fastq1 = None,
        fastq2 = None,
    ):
        self.alignment_file_path = alignment_file_path
        self.output_file_path = output_file_path
        self.protein_file_path = protein_file_path,
        self.nucleotide_file_path = nucleotide_file_path 
        self.alignment_list_path = alignment_list_path
        self.prefix = prefix
        self.fasta = fasta
        self.idmap = idmap
        self.reference = reference
        self.verbose = verbose
        self.entry = entry
        self.input_file = input_file
        self.output_file_format = output_file_format
        self.input_file_format = input_file_format
        self.output_file = output_file
        self.threshold = threshold
        self.ambiguous_character = ambiguous_character
        self.reverse = reverse
        self.translation_table = translation_table
        self.fastq = fastq
        self.minimum = minimum
        self.length = length
        self.output_file_1 = output_file_1
        self.output_file_2 = output_file_2
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def get_alignment_and_format(self):
        """
        automatic file type determination
        """
        try:
            return get_alignment_and_format_helper(self.alignment_file_path)
        except FileNotFoundError:
            print("Input corresponds to no such file or directory.")
            print("Please double check pathing and filenames")
            sys.exit()
    
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
        threshold = sum_contig_lens*.50
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
        threshold = sum_contig_lens*.90
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
        n50_threshold = sum_contig_lens*.50
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
        threshold = sum_contig_lens*.90
        n90 = 0

        for contig_len in contig_lens:
            n90 += contig_len
            if n90 >= threshold:
                return contig_len
                break
    
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
    
    def character_frequency(self) -> dict:
        """
        determine the frequency of characters in a
        FASTA file
        """
        # get contig lengths
        records = SeqIO.parse(self.fasta, "fasta")
        seqs = []
        for record in records:
            seqs.append(record.seq._data.upper())
        res = dict(Counter(''.join(seqs)))
        
        return res


