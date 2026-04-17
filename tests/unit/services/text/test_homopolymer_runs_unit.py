import json
from argparse import Namespace

from biokit.services.text.homopolymer_runs import HomopolymerRuns


def _args(fasta, per_base=False, format=None):
    return Namespace(fasta=fasta, per_base=per_base, format=format)


def test_homopolymer_runs_basic_tsv(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">seq1\nAAATTTTGCAAAAACGT\n>seq2\nGGGGCCCCAATAAA\n")
    HomopolymerRuns(_args(str(fa))).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0] == "id\tlength\tbase\tposition"
    # seq1: AAATTTTGCAAAAACGT — longest is 5 A's starting at position 10
    assert lines[1] == "seq1\t5\tA\t10"
    # seq2: GGGG tied with CCCC; first occurring is GGGG at position 1
    assert lines[2] == "seq2\t4\tG\t1"


def test_homopolymer_runs_is_case_insensitive(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\naaaaaccctttt\n")
    HomopolymerRuns(_args(str(fa))).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[1] == "s\t5\tA\t1"


def test_homopolymer_runs_json(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">a\nAATTTCCCC\n")
    HomopolymerRuns(_args(str(fa), format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert data == [{"id": "a", "length": 4, "base": "C", "position": 6}]


def test_homopolymer_runs_per_base(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nAATTTGCCCCGGG\n")
    HomopolymerRuns(_args(str(fa), per_base=True)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0] == "id\tA\tC\tG\tT"
    # A:2, C:4, G:3, T:3
    assert lines[1] == "s\t2\t4\t3\t3"


def test_homopolymer_runs_per_base_json(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nAAAAAGGGGTTTCC\n")
    HomopolymerRuns(_args(str(fa), per_base=True, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert data == [{"id": "s", "A": 5, "C": 2, "G": 4, "T": 3}]


def test_homopolymer_runs_single_character(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nA\n")
    HomopolymerRuns(_args(str(fa))).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[1] == "s\t1\tA\t1"


def test_homopolymer_runs_homopolymer_only(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nAAAAAAAA\n")
    HomopolymerRuns(_args(str(fa))).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[1] == "s\t8\tA\t1"


def test_homopolymer_runs_ties_return_first(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    # Both AAAA and TTTT are length 4; first occurrence wins
    fa.write_text(">s\nCCAAAAGGTTTT\n")
    HomopolymerRuns(_args(str(fa))).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[1] == "s\t4\tA\t3"


def test_homopolymer_runs_includes_ambiguous_bases(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">s\nATNNNNNAT\n")
    HomopolymerRuns(_args(str(fa))).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[1] == "s\t5\tN\t3"


def test_homopolymer_runs_per_base_ignores_ambiguous(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    # N run of 5 should not appear in per-base output (A/C/G/T only)
    fa.write_text(">s\nANNNNNACGT\n")
    HomopolymerRuns(_args(str(fa), per_base=True)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0] == "id\tA\tC\tG\tT"
    assert lines[1] == "s\t1\t1\t1\t1"
