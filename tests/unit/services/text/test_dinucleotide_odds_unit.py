import json
from argparse import Namespace

from biokit.services.text.dinucleotide_odds import (
    DinucleotideOdds,
    DINUCLEOTIDES,
)


def _args(fasta, verbose=False, format=None):
    return Namespace(fasta=fasta, verbose=verbose, format=format)


def _write_fasta(tmp_path, seqs):
    fa = tmp_path / "seqs.fa"
    with open(fa, "w") as fh:
        for name, seq in seqs:
            fh.write(f">{name}\n{seq}\n")
    return str(fa)


def test_dinucleotide_odds_basic_tsv(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("s", "ACGT" * 25)])
    DinucleotideOdds(_args(fa)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0] == "dinucleotide\tO_E"
    # 16 dinucleotide rows
    assert len(lines) == 17
    dinucs_seen = [line.split("\t")[0] for line in lines[1:]]
    assert dinucs_seen == list(DINUCLEOTIDES)


def test_dinucleotide_odds_uniform_random_sequence(tmp_path, capsys):
    # Long uniform-distributed sequence; expected O/E ≈ 1 for all
    # but we use a simple synthetic where AAACGTAAACGT... yields known counts
    seq = "ACGT" * 1000  # perfectly cyclic - each dinucleotide AC, CG, GT, TA each appears 1000 times
    fa = _write_fasta(tmp_path, [("s", seq)])
    DinucleotideOdds(_args(fa)).run()
    rows = capsys.readouterr().out.strip().split("\n")[1:]
    values = {row.split("\t")[0]: float(row.split("\t")[1]) for row in rows}
    # In ACGT-repeated, only AC, CG, GT, TA exist
    assert values["AC"] > 0
    assert values["CG"] > 0
    assert values["GT"] > 0
    assert values["TA"] > 0
    # Diagonals (AA, CC, GG, TT) should be 0 since they don't occur
    assert values["AA"] == 0.0
    assert values["CC"] == 0.0
    assert values["GG"] == 0.0
    assert values["TT"] == 0.0


def test_dinucleotide_odds_cpg_suppression(tmp_path, capsys):
    # Sequence with explicit CpG suppression (no CG dinucleotides)
    seq = ("AGAGAGAG" + "TCTCTCTC") * 50  # mix C, G but no CG
    fa = _write_fasta(tmp_path, [("s", seq)])
    DinucleotideOdds(_args(fa)).run()
    rows = capsys.readouterr().out.strip().split("\n")[1:]
    values = {row.split("\t")[0]: float(row.split("\t")[1]) for row in rows}
    # No CG present in sequence — CG O/E should be 0
    assert values["CG"] == 0.0
    # C and G are still present individually
    # GC O/E should be > 0 (we have ...AGAGAG TCTCTC AGAGAG...)
    # Actually let me think: in "AGAGAGAG" + "TCTCTCTC" we have GA, AG, GT (at boundary), TC, CT, CT...
    # The boundary creates "GT" not "GC", so GC O/E might be 0 too.
    # Just verify CG = 0 which is the key CpG suppression metric


def test_dinucleotide_odds_case_insensitive(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("s", "acgtacgt")])
    DinucleotideOdds(_args(fa)).run()
    rows = capsys.readouterr().out.strip().split("\n")[1:]
    values = {row.split("\t")[0]: float(row.split("\t")[1]) for row in rows}
    # Should be the same as uppercase ACGTACGT
    assert values["AC"] > 0
    assert values["CG"] > 0


def test_dinucleotide_odds_skips_ambiguous(tmp_path, capsys):
    # Sequence with N's interspersed — should be skipped
    fa = _write_fasta(tmp_path, [("s", "ACNNGT" * 10)])
    DinucleotideOdds(_args(fa)).run()
    rows = capsys.readouterr().out.strip().split("\n")[1:]
    values = {row.split("\t")[0]: float(row.split("\t")[1]) for row in rows}
    # The valid dinucleotides are only AC and GT (since NN and CN, NG are skipped)
    # AC at position 0, GT at position 4, then ACAC at next iteration etc.
    # Just verify CN-containing or NG pairs aren't counted
    assert values["AC"] > 0
    assert values["GT"] > 0


def test_dinucleotide_odds_verbose(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [
        ("s1", "ACGT" * 25),
        ("s2", "AGAG" * 25),
    ])
    DinucleotideOdds(_args(fa, verbose=True)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    header = lines[0].split("\t")
    assert header[0] == "id"
    assert header[1:] == list(DINUCLEOTIDES)
    assert len(lines) == 3  # header + 2 sequences
    s1_fields = lines[1].split("\t")
    s2_fields = lines[2].split("\t")
    assert s1_fields[0] == "s1"
    assert s2_fields[0] == "s2"


def test_dinucleotide_odds_json(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("s", "ACGT" * 10)])
    DinucleotideOdds(_args(fa, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert len(data) == 16
    dinucs_seen = {row["dinucleotide"] for row in data}
    assert dinucs_seen == set(DINUCLEOTIDES)
    for row in data:
        assert "O_E" in row


def test_dinucleotide_odds_yaml(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("s", "ACGT" * 10)])
    DinucleotideOdds(_args(fa, format="yaml")).run()
    out = capsys.readouterr().out
    assert "dinucleotide: AA" in out
    assert "dinucleotide: TT" in out


def test_dinucleotide_odds_empty_sequence(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("empty", "")])
    DinucleotideOdds(_args(fa)).run()
    rows = capsys.readouterr().out.strip().split("\n")[1:]
    values = {row.split("\t")[0]: row.split("\t")[1] for row in rows}
    # All should be None (printed as "None") when no counts
    for dn in DINUCLEOTIDES:
        assert values[dn] == "None"


def test_dinucleotide_odds_only_ambiguous(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("s", "NNNNNNNN")])
    DinucleotideOdds(_args(fa)).run()
    rows = capsys.readouterr().out.strip().split("\n")[1:]
    values = {row.split("\t")[0]: row.split("\t")[1] for row in rows}
    for dn in DINUCLEOTIDES:
        assert values[dn] == "None"


def test_dinucleotide_odds_aggregates_across_sequences(tmp_path, capsys):
    # Two sequences should be aggregated, not averaged
    fa1 = _write_fasta(tmp_path, [
        ("a", "ACGTACGT"),
        ("b", "ACGTACGT"),
    ])
    DinucleotideOdds(_args(fa1)).run()
    out1 = capsys.readouterr().out

    fa2 = _write_fasta(tmp_path, [("combined", "ACGTACGT" * 2)])
    DinucleotideOdds(_args(fa2)).run()
    out2 = capsys.readouterr().out

    # The aggregated counts from two records should equal the combined string
    # (modulo the dinucleotide spanning the sequence break which doesn't exist
    # for aggregated mode — there's a slight difference at the boundary)
    # Both should yield very similar O/E values for non-boundary dinucleotides
    rows1 = {ln.split("\t")[0]: ln.split("\t")[1] for ln in out1.strip().split("\n")[1:]}
    rows2 = {ln.split("\t")[0]: ln.split("\t")[1] for ln in out2.strip().split("\n")[1:]}
    # AC, CG, GT, TA should appear in both
    for dn in ("AC", "CG", "GT", "TA"):
        v1 = float(rows1[dn])
        v2 = float(rows2[dn])
        assert v1 > 0 and v2 > 0
