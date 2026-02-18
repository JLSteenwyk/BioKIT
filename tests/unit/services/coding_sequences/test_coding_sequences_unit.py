from argparse import Namespace

from biokit.services.coding_sequences.gc_content_first_position import GCContentFirstPosition
from biokit.services.coding_sequences.gc_content_second_position import GCContentSecondPosition
from biokit.services.coding_sequences.gc_content_third_position import GCContentThirdPosition
from biokit.services.coding_sequences.gene_wise_relative_synonymous_codon_usage import (
    GeneWiseRelativeSynonymousCodonUsage,
)


def test_gc_content_first_position_non_verbose_aggregates_all_records(tmp_path, capsys):
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">a\nGAA\n>b\nATA\n")
    args = Namespace(fasta=str(fasta), verbose=False)

    GCContentFirstPosition(args).run()
    assert capsys.readouterr().out.strip() == "0.5"


def test_gc_content_second_position_non_verbose_aggregates_all_records(tmp_path, capsys):
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">a\nAGA\n>b\nAAA\n")
    args = Namespace(fasta=str(fasta), verbose=False)

    GCContentSecondPosition(args).run()
    assert capsys.readouterr().out.strip() == "0.5"


def test_gc_content_third_position_non_verbose_aggregates_all_records(tmp_path, capsys):
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">a\nAAG\n>b\nAAA\n")
    args = Namespace(fasta=str(fasta), verbose=False)

    GCContentThirdPosition(args).run()
    assert capsys.readouterr().out.strip() == "0.5"


def test_gene_wise_rscu_single_codon_gene_does_not_crash(tmp_path, capsys):
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">gene1\nATG\n")
    args = Namespace(fasta=str(fasta), translation_table=None)

    GeneWiseRelativeSynonymousCodonUsage(args).run()
    assert capsys.readouterr().out.strip() == "gene1\t1.0\t1.0\t0.0"
