import json
from argparse import Namespace

from biokit.services.alignment.transition_transversion_ratio import (
    TransitionTransversionRatio,
)


def _args(fasta, verbose=False, format=None):
    return Namespace(fasta=fasta, verbose=verbose, format=format)


def _write_aln(tmp_path, seqs):
    fa = tmp_path / "aln.fa"
    with open(fa, "w") as fh:
        for ident, s in seqs:
            fh.write(f">{ident}\n{s}\n")
    return str(fa)


def test_ti_tv_pure_transitions(tmp_path, capsys):
    # 2 sequences differing only by A<->G transitions
    aln = _write_aln(tmp_path, [
        ("s1", "AAAA"),
        ("s2", "GGGG"),
    ])
    TransitionTransversionRatio(_args(aln)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    assert fields[0] == "4"   # 4 transitions
    assert fields[1] == "0"   # 0 transversions
    assert fields[2] == "None"  # ratio undefined (tv == 0)


def test_ti_tv_pure_transversions(tmp_path, capsys):
    # 2 sequences differing only by A<->C transversions
    aln = _write_aln(tmp_path, [
        ("s1", "AAAA"),
        ("s2", "CCCC"),
    ])
    TransitionTransversionRatio(_args(aln)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    assert fields[0] == "0"
    assert fields[1] == "4"
    assert float(fields[2]) == 0.0


def test_ti_tv_mixed_ratio(tmp_path, capsys):
    # Hand-built alignment with known counts
    # col1: A,A → constant
    # col2: A,G → 1 Ti
    # col3: A,C → 1 Tv
    # col4: C,T → 1 Ti (pyrimidine)
    # col5: G,T → 1 Tv (purine-pyrimidine)
    aln = _write_aln(tmp_path, [
        ("s1", "AAACG"),
        ("s2", "AGCTT"),
    ])
    TransitionTransversionRatio(_args(aln)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    assert fields[0] == "2"   # transitions: A↔G, C↔T
    assert fields[1] == "2"   # transversions: A↔C, G↔T
    assert float(fields[2]) == 1.0


def test_ti_tv_skips_gap_columns(tmp_path, capsys):
    # Column 2 contains a gap → skipped entirely
    aln = _write_aln(tmp_path, [
        ("s1", "AAG"),
        ("s2", "G-T"),
    ])
    TransitionTransversionRatio(_args(aln)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    # col1: A↔G = 1 Ti
    # col2: contains gap → skipped
    # col3: G↔T = 1 Tv
    assert fields[0] == "1"
    assert fields[1] == "1"
    assert float(fields[2]) == 1.0


def test_ti_tv_question_mark_treated_as_gap(tmp_path, capsys):
    aln = _write_aln(tmp_path, [
        ("s1", "AGC"),
        ("s2", "AG?"),
    ])
    TransitionTransversionRatio(_args(aln)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    # col1: constant, col2: constant, col3: '?' → gap
    assert fields[0] == "0"
    assert fields[1] == "0"
    assert fields[2] == "None"


def test_ti_tv_verbose_classification(tmp_path, capsys):
    aln = _write_aln(tmp_path, [
        ("s1", "AAGCG"),
        ("s2", "AGCC-"),
    ])
    TransitionTransversionRatio(_args(aln, verbose=True)).run()
    out = capsys.readouterr().out.strip().split("\n")
    # col1: A,A → constant
    # col2: A,G → transition
    # col3: G,C → transversion
    # col4: C,C → constant
    # col5: G,- → gap
    assert out[0] == "1\tconstant"
    assert out[1] == "2\ttransition"
    assert out[2] == "3\ttransversion"
    assert out[3] == "4\tconstant"
    assert out[4] == "5\tgap"


def test_ti_tv_multi_sequence_pairwise_count(tmp_path, capsys):
    # 4 sequences, 1 column: A,A,G,G
    # Pairs (A,G) appear 4 times (1-3, 1-4, 2-3, 2-4) → 4 transitions
    aln = _write_aln(tmp_path, [
        ("s1", "A"),
        ("s2", "A"),
        ("s3", "G"),
        ("s4", "G"),
    ])
    TransitionTransversionRatio(_args(aln)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    assert fields[0] == "4"
    assert fields[1] == "0"


def test_ti_tv_verbose_marks_mixed_column_as_transversion(tmp_path, capsys):
    # Column with A, G, C — has both a transition (A<->G) and
    # transversions (A<->C, G<->C). Label should be 'transversion'.
    aln = _write_aln(tmp_path, [
        ("s1", "A"),
        ("s2", "G"),
        ("s3", "C"),
    ])
    TransitionTransversionRatio(_args(aln, verbose=True)).run()
    assert capsys.readouterr().out.strip() == "1\ttransversion"


def test_ti_tv_json_summary(tmp_path, capsys):
    aln = _write_aln(tmp_path, [
        ("s1", "AAAA"),
        ("s2", "GCGT"),
    ])
    TransitionTransversionRatio(_args(aln, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert set(data.keys()) == {"transitions", "transversions", "ratio"}
    # col1: A↔G Ti, col2: A↔C Tv, col3: A↔G Ti, col4: A↔T Tv
    assert data["transitions"] == 2
    assert data["transversions"] == 2
    assert data["ratio"] == 1.0


def test_ti_tv_json_verbose(tmp_path, capsys):
    aln = _write_aln(tmp_path, [
        ("s1", "AGC"),
        ("s2", "AGT"),
    ])
    TransitionTransversionRatio(_args(aln, verbose=True, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert len(data) == 3
    assert data[0] == {"position": 1, "classification": "constant"}
    assert data[1] == {"position": 2, "classification": "constant"}
    assert data[2] == {"position": 3, "classification": "transition"}


def test_ti_tv_case_insensitive(tmp_path, capsys):
    aln = _write_aln(tmp_path, [
        ("s1", "aaaa"),
        ("s2", "ggcc"),
    ])
    TransitionTransversionRatio(_args(aln)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    # cols 1,2: a↔g (2 Ti); cols 3,4: a↔c (2 Tv)
    assert fields[0] == "2"
    assert fields[1] == "2"
    assert float(fields[2]) == 1.0
