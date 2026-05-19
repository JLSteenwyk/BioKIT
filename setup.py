from os import path
from pathlib import Path
import re
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))
version_text = Path(here, "biokit", "version.py").read_text(encoding="utf-8")
version_match = re.search(r'__version__\s*=\s*"([^"]+)"', version_text)
if version_match is None:
    raise RuntimeError("Unable to find __version__ in biokit/version.py")
__version__ = version_match.group(1)

with open(path.join(here, "README.md")) as f:
    long_description = f.read()

CLASSIFIERS = [
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Topic :: Scientific/Engineering",
]

REQUIRES = ["biopython>=1.82", "matplotlib>=3.6.0", "numpy>=1.24.0", "cython"]

setup(
    name="jlsteenwyk-biokit",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/biokit",
    packages=find_packages(),
    classifiers=CLASSIFIERS,
    entry_points={
        "console_scripts": [
            "biokit = biokit.biokit:main",
            "bk_alignment_length = biokit.biokit:alignment_length",  # Alignment-based functions
            "bk_aln_len = biokit.biokit:alignment_length",
            "bk_alignment_recoding = biokit.biokit:alignment_recoding",
            "bk_aln_recoding = biokit.biokit:alignment_recoding",
            "bk_recode = biokit.biokit:alignment_recoding",
            "bk_alignment_summary = biokit.biokit:alignment_summary",
            "bk_aln_summary = biokit.biokit:alignment_summary",
            "bk_consensus_sequence = biokit.biokit:consensus_sequence",
            "bk_con_seq = biokit.biokit:consensus_sequence",
            "bk_constant_sites = biokit.biokit:constant_sites",
            "bk_con_sites = biokit.biokit:constant_sites",
            "bk_parsimony_informative_sites = biokit.biokit:parsimony_informative_sites",
            "bk_pi_sites = biokit.biokit:parsimony_informative_sites",
            "bk_pis = biokit.biokit:parsimony_informative_sites",
            "bk_position_specific_score_matrix = biokit.biokit:position_specific_score_matrix",
            "bk_pssm = biokit.biokit:position_specific_score_matrix",
            "bk_transition_transversion_ratio = biokit.biokit:transition_transversion_ratio",
            "bk_ti_tv = biokit.biokit:transition_transversion_ratio",
            "bk_variable_sites = biokit.biokit:variable_sites",
            "bk_var_sites = biokit.biokit:variable_sites",
            "bk_vs = biokit.biokit:variable_sites",
            "bk_gc_content_first_position = biokit.biokit:gc_content_first_position",  # Coding sequences-based functions
            "bk_gc1 = biokit.biokit:gc_content_first_position",
            "bk_gc_content_second_position = biokit.biokit:gc_content_second_position",
            "bk_gc2 = biokit.biokit:gc_content_second_position",
            "bk_gc_content_third_position = biokit.biokit:gc_content_third_position",
            "bk_gc3 = biokit.biokit:gc_content_third_position",
            "bk_gc_content_four_fold_degenerate_sites = biokit.biokit:gc_content_four_fold_degenerate_sites",
            "bk_gc4 = biokit.biokit:gc_content_four_fold_degenerate_sites",
            "bk_neutrality_plot = biokit.biokit:neutrality_plot",
            "bk_effective_number_of_codons = biokit.biokit:effective_number_of_codons",
            "bk_enc = biokit.biokit:effective_number_of_codons",
            "bk_codon_adaptation_index = biokit.biokit:codon_adaptation_index",
            "bk_cai = biokit.biokit:codon_adaptation_index",
            "bk_gene_wise_relative_synonymous_codon_usage = biokit.biokit:gene_wise_relative_synonymous_codon_usage",
            "bk_gene_wise_rscu = biokit.biokit:gene_wise_relative_synonymous_codon_usage",
            "bk_gw_rscu = biokit.biokit:gene_wise_relative_synonymous_codon_usage",
            "bk_grscu = biokit.biokit:gene_wise_relative_synonymous_codon_usage",
            "bk_relative_synonymous_codon_usage = biokit.biokit:relative_synonymous_codon_usage",
            "bk_rscu = biokit.biokit:relative_synonymous_codon_usage",
            "bk_translate_sequence = biokit.biokit:translate_sequence",
            "bk_translate_seq = biokit.biokit:translate_sequence",
            "bk_trans_seq = biokit.biokit:translate_sequence",
            "bk_fastq_read_lengths = biokit.biokit:fastq_read_lengths",  # FASTQ-based functions
            "bk_fastq_read_lens = biokit.biokit:fastq_read_lengths",
            "bk_fastq_to_fasta = biokit.biokit:fastq_to_fasta",
            "bk_fq2fa = biokit.biokit:fastq_to_fasta",
            "bk_subset_pe_fastq_reads = biokit.biokit:subset_pe_fastq_reads",
            "bk_subset_pe_fastq = biokit.biokit:subset_pe_fastq_reads",
            "bk_subset_se_fastq_reads = biokit.biokit:subset_se_fastq_reads",
            "bk_subset_se_fastq = biokit.biokit:subset_se_fastq_reads",
            "bk_trim_pe_adapters_fastq = biokit.biokit:trim_pe_adapters_fastq",
            "bk_trim_pe_adapters_fastq_reads = biokit.biokit:trim_pe_adapters_fastq",
            "bk_trim_pe_fastq = biokit.biokit:trim_pe_fastq",
            "bk_trim_pe_fastq_reads = biokit.biokit:trim_pe_fastq",
            "bk_trim_se_adapters_fastq = biokit.biokit:trim_se_adapters_fastq",
            "bk_trim_se_adapters_fastq_reads = biokit.biokit:trim_se_adapters_fastq",
            "bk_trim_se_fastq = biokit.biokit:trim_se_fastq",
            "bk_trim_se_fastq_reads = biokit.biokit:trim_se_fastq",
            "bk_assembly_curve = biokit.biokit:assembly_curve",  # genomes-based functions
            "bk_asm_curve = biokit.biokit:assembly_curve",
            "bk_gc_content = biokit.biokit:gc_content",
            "bk_gc = biokit.biokit:gc_content",
            "bk_genome_assembly_metrics = biokit.biokit:genome_assembly_metrics",
            "bk_assembly_metrics = biokit.biokit:genome_assembly_metrics",
            "bk_l50 = biokit.biokit:l50",
            "bk_l90 = biokit.biokit:l90",
            "bk_longest_scaffold = biokit.biokit:longest_scaffold",
            "bk_longest_scaff = biokit.biokit:longest_scaffold",
            "bk_longest_contig = biokit.biokit:longest_scaffold",
            "bk_longest_cont = biokit.biokit:longest_scaffold",
            "bk_n50 = biokit.biokit:n50",
            "bk_n90 = biokit.biokit:n90",
            "bk_number_of_large_scaffolds = biokit.biokit:number_of_large_scaffolds",
            "bk_num_of_lrg_scaffolds = biokit.biokit:number_of_large_scaffolds",
            "bk_number_of_large_contigs = biokit.biokit:number_of_large_scaffolds",
            "bk_num_of_lrg_cont = biokit.biokit:number_of_large_scaffolds",
            "bk_number_of_scaffolds = biokit.biokit:number_of_scaffolds",
            "bk_num_of_scaffolds = biokit.biokit:number_of_scaffolds",
            "bk_number_of_contigs = biokit.biokit:number_of_scaffolds",
            "bk_num_of_cont = biokit.biokit:number_of_scaffolds",
            "bk_sum_of_scaffold_lengths = biokit.biokit:sum_of_scaffold_lengths",
            "bk_sum_of_contig_lengths = biokit.biokit:sum_of_scaffold_lengths",
            "bk_character_frequency = biokit.biokit:character_frequency",  # text-based functions
            "bk_char_freq = biokit.biokit:character_frequency",
            "bk_dinucleotide_odds = biokit.biokit:dinucleotide_odds",
            "bk_dno = biokit.biokit:dinucleotide_odds",
            "bk_faidx = biokit.biokit:faidx",
            "bk_get_entry = biokit.biokit:faidx",
            "bk_ge = biokit.biokit:faidx",
            "bk_file_format_converter = biokit.biokit:file_format_converter",
            "bk_format_converter = biokit.biokit:file_format_converter",
            "bk_ffc = biokit.biokit:file_format_converter",
            "bk_genbank_to_fasta = biokit.biokit:genbank_to_fasta",
            "bk_gb2fa = biokit.biokit:genbank_to_fasta",
            "bk_homopolymer_runs = biokit.biokit:homopolymer_runs",
            "bk_homopolymer = biokit.biokit:homopolymer_runs",
            "bk_kmer_frequency = biokit.biokit:kmer_frequency",
            "bk_kmer_freq = biokit.biokit:kmer_frequency",
            "bk_find_orfs = biokit.biokit:find_orfs",
            "bk_orfs = biokit.biokit:find_orfs",
            "bk_multiple_line_to_single_line_fasta = biokit.biokit:multiple_line_to_single_line_fasta",
            "bk_ml2sl = biokit.biokit:multiple_line_to_single_line_fasta",
            "bk_remove_short_sequences = biokit.biokit:remove_short_sequences",
            "bk_remove_short_seqs = biokit.biokit:remove_short_sequences",
            "bk_fasta_deduplication = biokit.biokit:fasta_deduplication",
            "bk_dedup = biokit.biokit:fasta_deduplication",
            "bk_melting_temperature = biokit.biokit:melting_temperature",
            "bk_tm = biokit.biokit:melting_temperature",
            "bk_protein_charge = biokit.biokit:protein_charge",
            "bk_prot_charge = biokit.biokit:protein_charge",
            "bk_protein_properties = biokit.biokit:protein_properties",
            "bk_prot_prop = biokit.biokit:protein_properties",
            "bk_remove_fasta_entry = biokit.biokit:remove_fasta_entry",
            "bk_rename_fasta_entries = biokit.biokit:rename_fasta_entries",
            "bk_rename_fasta = biokit.biokit:rename_fasta_entries",
            "bk_restriction_sites = biokit.biokit:restriction_sites",
            "bk_re_sites = biokit.biokit:restriction_sites",
            "bk_sample_sequences = biokit.biokit:sample_sequences",
            "bk_sample_seqs = biokit.biokit:sample_sequences",
            "bk_shuffle_sequences = biokit.biokit:shuffle_sequences",
            "bk_shuffle_seqs = biokit.biokit:shuffle_sequences",
            "bk_reorder_by_sequence_length = biokit.biokit:reorder_by_sequence_length",
            "bk_reorder_by_seq_len = biokit.biokit:reorder_by_sequence_length",
            "bk_sequence_complement = biokit.biokit:sequence_complement",
            "bk_seq_comp = biokit.biokit:sequence_complement",
            "bk_sequence_length = biokit.biokit:sequence_length",
            "bk_seq_len = biokit.biokit:sequence_length",
            "bk_single_line_to_multiple_line_fasta = biokit.biokit:single_line_to_multiple_line_fasta",
            "bk_sl2ml = biokit.biokit:single_line_to_multiple_line_fasta",
        ]
    },
    version=__version__,
    python_requires=">=3.11",
    include_package_data=True,
    install_requires=REQUIRES,
)

# push new version to pypi
# rm -rf dist
# python setup.py sdist bdist_wheel --universal
# twine upload dist/* -r pypi
# then push to anaconda
