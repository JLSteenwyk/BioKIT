from argparse import Namespace

from biokit.services.coding_sequences.gc_content_four_fold_degenerate_sites import (
    GCContentFourFoldDegenerateSites,
)


def test_gc4_non_verbose_aggregates(tmp_path, capsys):
    # GCU=A (four-fold), GCG=A (four-fold) -> third positions: U, G
    # In DNA: GCT, GCG -> third positions: T, G -> GC = 1/2 = 0.5
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">a\nGCTGCG\n")
    args = Namespace(fasta=str(fasta), verbose=False, translation_table=None)

    GCContentFourFoldDegenerateSites(args).run()
    assert capsys.readouterr().out.strip() == "0.5"


def test_gc4_verbose_per_sequence(tmp_path, capsys):
    # seq a: GCT (A, 4-fold, 3rd=T), GCG (A, 4-fold, 3rd=G) -> GC = 0.5
    # seq b: GCC (A, 4-fold, 3rd=C), GCC (A, 4-fold, 3rd=C) -> GC = 1.0
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">a\nGCTGCG\n>b\nGCCGCC\n")
    args = Namespace(fasta=str(fasta), verbose=True, translation_table=None)

    GCContentFourFoldDegenerateSites(args).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 2
    assert lines[0] == "a\t0.5"
    assert lines[1] == "b\t1.0"


def test_gc4_skips_non_four_fold_codons(tmp_path, capsys):
    # AAA=K (not four-fold: AAA=K, AAG=K, AAU=N, AAC=N)
    # GCT=A (four-fold, 3rd=T)
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">a\nAAAGCT\n")
    args = Namespace(fasta=str(fasta), verbose=False, translation_table=None)

    GCContentFourFoldDegenerateSites(args).run()
    # Only one four-fold site (GCT), third pos = T -> GC = 0.0
    assert capsys.readouterr().out.strip() == "0.0"


def test_gc4_multiple_sequences_aggregate(tmp_path, capsys):
    # seq a: GCG (4-fold, 3rd=G), seq b: GCA (4-fold, 3rd=A)
    # Aggregated third positions: G, A -> GC = 1/2 = 0.5
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">a\nGCG\n>b\nGCA\n")
    args = Namespace(fasta=str(fasta), verbose=False, translation_table=None)

    GCContentFourFoldDegenerateSites(args).run()
    assert capsys.readouterr().out.strip() == "0.5"


def test_gc4_all_four_fold_families(tmp_path, capsys):
    # One codon from each 4-fold degenerate family (standard code):
    # CTG=L, GTG=V, TCG=S, CCG=P, ACG=T, GCG=A, CGG=R, GGG=G
    # All third positions are G -> GC = 1.0
    fasta = tmp_path / "cds.fa"
    fasta.write_text(">a\nCTGGTGTCGCCGACGGCGCGGGGG\n")
    args = Namespace(fasta=str(fasta), verbose=False, translation_table=None)

    GCContentFourFoldDegenerateSites(args).run()
    assert capsys.readouterr().out.strip() == "1.0"
