from collections import Counter
from os import path
import sys

from Bio import SeqIO

from ..base import BaseService
from ...helpers.files import get_alignment_and_format as get_alignment_and_format_helper

here = path.dirname(__file__)

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
        percent = None,
        seed = None,
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
        self.percent = percent
        self.seed = seed

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

    def read_translation_table(self) -> dict:
        """
        return translation table with codons as keys and amino acids as values
        """

        trans_table = dict()
        
        if self.translation_table is None or self.translation_table == '1':
            pathing = path.join(here, "../../tables/standard_genetic_code.txt")
        elif self.translation_table == '2':
            pathing = path.join(here, "../../tables/vertebrate_mitochondrial_code.txt")
        elif self.translation_table == '3':
            pathing = path.join(here, "../../tables/yeast_mitochondrial_code.txt")
        elif self.translation_table == '4':
            pathing = path.join(here, "../../tables/mold_protozoan_and_coelenterate_mitochondrial_code_and_the_mycoplasma_spiroplasma.txt")
        elif self.translation_table == '5':
            pathing = path.join(here, "../../tables/invertebrate_mitochondrial_code.txt")
        elif self.translation_table == '6':
            pathing = path.join(here, "../../tables/ciliate_dasycladacean_and_hexamita_nuclear_code.txt")
        elif self.translation_table == '9':
            pathing = path.join(here, "../../tables/echinoderm_and_flatworm_mitochondrial_code.txt")
        elif self.translation_table == '10':
            pathing = path.join(here, "../../tables/euplotid_nuclear_code.txt")
        elif self.translation_table == '11':
            pathing = path.join(here, "../../tables/bacterial_archaeal_and_plant_plastid_code.txt")
        elif self.translation_table == '12':
            pathing = path.join(here, "../../tables/alternative_yeast_nuclear_code.txt")
        elif self.translation_table == '13':
            pathing = path.join(here, "../../tables/ascidian_mitochondrial_code.txt")
        elif self.translation_table == '14':
            pathing = path.join(here, "../../tables/alternative_flatworm_mitochondrial_code.txt")
        elif self.translation_table == '16':
            pathing = path.join(here, "../../tables/chlorophycean_mitochondrial_code.txt")
        elif self.translation_table == '21':
            pathing = path.join(here, "../../tables/trematode_mitochondrial_code.txt")
        elif self.translation_table == '22':
            pathing = path.join(here, "../../tables/scenedesmus_obliquus_mitochondrial_code.txt")
        elif self.translation_table == '23':
            pathing = path.join(here, "../../tables/thraustochytrium_mitochondrial_code.txt")
        elif self.translation_table == '24':
            pathing = path.join(here, "../../tables/rhabdopleuridae_mitochondrial_code.txt")
        elif self.translation_table == '25':
            pathing = path.join(here, "../../tables/candidate_division_sr1_and_gracilibacteria_code.txt")
        elif self.translation_table == '26':
            pathing = path.join(here, "../../tables/pachysolen_tannophilus_nuclear_code.txt")
        elif self.translation_table == '27':
            pathing = path.join(here, "../../tables/karyorelict_nuclear_code.txt")
        elif self.translation_table == '28':
            pathing = path.join(here, "../../tables/condylostoma_nuclear_code.txt")
        elif self.translation_table == '29':
            pathing = path.join(here, "../../tables/mesodinium_nuclear_code.txt")
        elif self.translation_table == '30':
            pathing = path.join(here, "../../tables/peritrich_nuclear_code.txt")
        elif self.translation_table == '31':
            pathing = path.join(here, "../../tables/blastocrithidia_nuclear_code.txt")
        elif self.translation_table == '33':
            pathing = path.join(here, "../../tables/cephalodiscidae_mitochondrial_UAA_tyr_code.txt")
        elif self.translation_table == '50':
            pathing = path.join(here, "../../tables/CUG_ala_code.txt")
        # case handling for a custom translation table
        else:
            pathing = str(pathing)

        with open(pathing) as code:
            for line in code:
                line=line.split()
                trans_table[line[0]] = line[1]
        return trans_table


