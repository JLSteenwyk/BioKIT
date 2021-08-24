import os
import pytest
from mock import patch, call # noqa
from pathlib import Path

from biokit.biokit import Biokit # noqa

here = Path(__file__)


@pytest.mark.slow
@pytest.mark.integration
class TestCLIRunners(object):
    def test_bk_alignment_summary(self):
        cmd = "bk_alignment_summary -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_aln_summary(self):
        cmd = "bk_aln_summary -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_consensus_sequence(self):
        cmd = "bk_consensus_sequence -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_con_seq(self):
        cmd = "bk_con_seq -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_position_specific_score_matrix(self):
        cmd = "bk_position_specific_score_matrix -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_pssm(self):
        cmd = "bk_pssm -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_gc_content_first_position(self):
        cmd = "bk_gc_content_first_position -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_gc1(self):
        cmd = "bk_gc1 -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_gc_content_second_position(self):
        cmd = "bk_gc_content_second_position -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_gc2(self):
        cmd = "bk_gc2 -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_gc_content_third_position(self):
        cmd = "bk_gc_content_third_position -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_gc3(self):
        cmd = "bk_gc3 -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_relative_synonymous_codon_usage(self):
        cmd = "bk_relative_synonymous_codon_usage -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_rscu(self):
        cmd = "bk_rscu -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_translate_sequence(self):
        cmd = "bk_translate_sequence -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_translate_seq(self):
        cmd = "bk_translate_seq -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_trans_seq(self):
        cmd = "bk_trans_seq -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_fastq_read_lengths(self):
        cmd = "bk_fastq_read_lengths -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_fastq_read_lens(self):
        cmd = "bk_fastq_read_lens -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_subset_pe_fastq_reads(self):
        cmd = "bk_subset_pe_fastq_reads -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_subset_pe_fastq(self):
        cmd = "bk_subset_pe_fastq -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_trim_pe_fastq_reads(self):
        cmd = "bk_trim_pe_fastq_reads -h"
        exit_status = os.system(cmd)
        assert exit_status == 32512

    def test_bk_trim_se_fastq_reads(self):
        cmd = "bk_trim_se_fastq_reads -h"
        exit_status = os.system(cmd)
        assert exit_status == 32512

    def test_bk_subset_se_fastq_reads(self):
        cmd = "bk_subset_se_fastq_reads -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_subset_se_fastq(self):
        cmd = "bk_subset_se_fastq -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_gc_content(self):
        cmd = "bk_gc_content -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_gc(self):
        cmd = "bk_gc -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_genome_assembly_metrics(self):
        cmd = "bk_genome_assembly_metrics -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_assembly_metrics(self):
        cmd = "bk_assembly_metrics -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_longest_scaffold(self):
        cmd = "bk_longest_scaffold -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_longest_scaff(self):
        cmd = "bk_longest_scaff -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_longest_contig(self):
        cmd = "bk_longest_contig -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_longest_cont(self):
        cmd = "bk_longest_cont -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_number_of_large_scaffolds(self):
        cmd = "bk_number_of_large_scaffolds -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_num_of_lrg_scaffolds(self):
        cmd = "bk_num_of_lrg_scaffolds -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_number_of_large_contigs(self):
        cmd = "bk_number_of_large_contigs -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_num_of_lrg_cont(self):
        cmd = "bk_num_of_lrg_cont -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_number_of_scaffolds(self):
        cmd = "bk_number_of_scaffolds -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_num_of_scaffolds(self):
        cmd = "bk_num_of_scaffolds -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_number_of_contigs(self):
        cmd = "bk_number_of_contigs -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_num_of_cont(self):
        cmd = "bk_num_of_cont -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_sum_of_scaffold_lengths(self):
        cmd = "bk_sum_of_scaffold_lengths -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_sum_of_contig_lengths(self):
        cmd = "bk_sum_of_contig_lengths -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_character_frequency(self):
        cmd = "bk_character_frequency -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_char_freq(self):
        cmd = "bk_char_freq -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_faidx(self):
        cmd = "bk_faidx -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    def test_bk_get_entry(self):
        cmd = "bk_get_entry -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    def test_bk_ge(self):
        cmd = "bk_ge -h"
        exit_status = os.system(cmd)
        assert exit_status == 0

    def test_bk_file_format_converter(self):
        cmd = "bk_file_format_converter -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_format_converter(self):
        cmd = "bk_format_converter -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_ffc(self):
        cmd = "bk_ffc -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_multiple_line_to_single_line_fasta(self):
        cmd = "bk_multiple_line_to_single_line_fasta -h"
        exit_status = os.system(cmd)
        assert exit_status == 32512

    def test_bk_ml2sl(self):
        cmd = "bk_ml2sl -h"
        exit_status = os.system(cmd)
        assert exit_status == 32512

    def test_bk_rename_fasta_entries(self):
        cmd = "bk_rename_fasta_entries -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_rename_fasta(self):
        cmd = "bk_rename_fasta -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_reorder_by_sequence_length(self):
        cmd = "bk_reorder_by_sequence_length -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_reorder_by_seq_len(self):
        cmd = "bk_reorder_by_seq_len -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_sequence_complement(self):
        cmd = "bk_sequence_complement -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_seq_comp(self):
        cmd = "bk_seq_comp -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_sequence_length(self):
        cmd = "bk_sequence_length -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_seq_len(self):
        cmd = "bk_seq_len -h"
        exit_status = os.system(cmd)
        assert exit_status == 256

    def test_bk_single_line_to_multiple_line_fasta(self):
        cmd = "bk_single_line_to_multiple_line_fasta -h"
        exit_status = os.system(cmd)
        assert exit_status == 32512

    def test_bk_sl2ml(self):
        cmd = "bk_sl2ml -h"
        exit_status = os.system(cmd)
        assert exit_status == 32512
