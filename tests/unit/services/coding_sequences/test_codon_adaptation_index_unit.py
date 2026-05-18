import json
from argparse import Namespace

import pytest

from biokit.services.coding_sequences.codon_adaptation_index import (
    CodonAdaptationIndex,
)


def _args(fasta, reference, translation_table=None, verbose=False, format=None):
    return Namespace(
        fasta=fasta,
        reference=reference,
        translation_table=translation_table,
        verbose=verbose,
        format=format,
    )


def _write_fasta(tmp_path, name, seqs):
    fa = tmp_path / name
    with open(fa, "w") as fh:
        for ident, s in seqs:
            fh.write(f">{ident}\n{s}\n")
    return str(fa)


# A reference set where only CTG (Leu) and GCT (Ala) are ever used.
# Met (ATG) and stop (TAA) are not informative for CAI.
REFERENCE_BLOCK = ("ATG" + "CTG" * 20 + "GCT" * 20 + "TAA")


def test_cai_matched_query_is_one(tmp_path, capsys):
    ref = _write_fasta(tmp_path, "ref.fa", [
        ("r1", REFERENCE_BLOCK),
        ("r2", REFERENCE_BLOCK),
    ])
    query = _write_fasta(tmp_path, "query.fa", [
        ("matched", "ATG" + "CTG" * 20 + "GCT" * 20 + "TAA"),
    ])
    CodonAdaptationIndex(_args(query, ref)).run()
    line = capsys.readouterr().out.strip()
    gene, cai = line.split("\t")
    assert gene == "matched"
    assert float(cai) == 1.0


def test_cai_anti_matched_query_is_low(tmp_path, capsys):
    ref = _write_fasta(tmp_path, "ref.fa", [
        ("r1", REFERENCE_BLOCK),
    ])
    # TTA = Leu (rare), GCG = Ala (rare). Both should have w = pseudocount/max.
    query = _write_fasta(tmp_path, "query.fa", [
        ("rare", "ATG" + "TTA" * 20 + "GCG" * 20 + "TAA"),
    ])
    CodonAdaptationIndex(_args(query, ref)).run()
    _, cai = capsys.readouterr().out.strip().split("\t")
    assert float(cai) < 0.05


def test_cai_skips_sequences_not_divisible_by_3(tmp_path, capsys):
    ref = _write_fasta(tmp_path, "ref.fa", [("r1", REFERENCE_BLOCK)])
    query = _write_fasta(tmp_path, "query.fa", [
        ("good", "ATG" + "CTG" * 5 + "TAA"),  # 21 nt
        ("bad", "ATG" + "CTG" * 5 + "TAA" + "A"),  # 22 nt
    ])
    CodonAdaptationIndex(_args(query, ref)).run()
    out = capsys.readouterr().out.strip()
    lines = out.split("\n")
    assert len(lines) == 1
    assert lines[0].startswith("good")


def test_cai_verbose_includes_n_codons(tmp_path, capsys):
    ref = _write_fasta(tmp_path, "ref.fa", [("r1", REFERENCE_BLOCK)])
    query = _write_fasta(tmp_path, "query.fa", [
        ("g", "ATG" + "CTG" * 5 + "GCT" * 5 + "TAA"),
    ])
    CodonAdaptationIndex(_args(query, ref, verbose=True)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    assert fields[0] == "g"
    assert float(fields[1]) == 1.0
    # 5 Leu codons + 5 Ala codons = 10 scored codons. Met and stop excluded.
    assert int(fields[2]) == 10


def test_cai_intermediate_value(tmp_path, capsys):
    # Reference: CTG and GCT dominate but TTA and GCG also appear sometimes
    # so their w-values are between 0 and 1.
    ref_seq = (
        "ATG"
        + "CTG" * 30  # major Leu
        + "TTA" * 10  # minor Leu
        + "GCT" * 30  # major Ala
        + "GCG" * 10  # minor Ala
        + "TAA"
    )
    ref = _write_fasta(tmp_path, "ref.fa", [("r1", ref_seq)])
    query = _write_fasta(tmp_path, "query.fa", [
        ("major", "ATG" + "CTG" * 10 + "GCT" * 10 + "TAA"),
        ("mixed", "ATG" + "CTG" * 5 + "TTA" * 5 + "GCT" * 5 + "GCG" * 5 + "TAA"),
        ("minor", "ATG" + "TTA" * 10 + "GCG" * 10 + "TAA"),
    ])
    CodonAdaptationIndex(_args(query, ref)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    values = {ln.split("\t")[0]: float(ln.split("\t")[1]) for ln in lines}
    assert values["major"] == 1.0
    assert 0.0 < values["minor"] < values["mixed"] < values["major"]


def test_cai_json(tmp_path, capsys):
    ref = _write_fasta(tmp_path, "ref.fa", [("r1", REFERENCE_BLOCK)])
    query = _write_fasta(tmp_path, "query.fa", [
        ("g1", "ATG" + "CTG" * 5 + "TAA"),
        ("g2", "ATG" + "GCT" * 5 + "TAA"),
    ])
    CodonAdaptationIndex(_args(query, ref, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert len(data) == 2
    for row in data:
        assert set(row.keys()) == {"gene_id", "CAI"}


def test_cai_json_verbose(tmp_path, capsys):
    ref = _write_fasta(tmp_path, "ref.fa", [("r1", REFERENCE_BLOCK)])
    query = _write_fasta(tmp_path, "query.fa", [
        ("g1", "ATG" + "CTG" * 5 + "TAA"),
    ])
    CodonAdaptationIndex(_args(query, ref, verbose=True, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert set(data[0].keys()) == {"gene_id", "CAI", "n_codons"}


def test_cai_requires_reference(tmp_path):
    query = _write_fasta(tmp_path, "query.fa", [("g", "ATGTAA")])
    with pytest.raises(ValueError):
        CodonAdaptationIndex(_args(query, "")).run()


def test_cai_custom_translation_table(tmp_path, capsys):
    # bacterial code (11) — should also yield CAI=1 when matched
    ref = _write_fasta(tmp_path, "ref.fa", [("r1", REFERENCE_BLOCK)])
    query = _write_fasta(tmp_path, "query.fa", [
        ("matched", "ATG" + "CTG" * 5 + "GCT" * 5 + "TAA"),
    ])
    CodonAdaptationIndex(_args(query, ref, translation_table="11")).run()
    _, cai = capsys.readouterr().out.strip().split("\t")
    assert float(cai) == 1.0


def test_cai_no_scored_codons_returns_none(tmp_path, capsys):
    # Query made up entirely of single-codon AAs (Met) plus a stop.
    # No synonymous codons → no scored codons → CAI undefined.
    ref = _write_fasta(tmp_path, "ref.fa", [("r1", REFERENCE_BLOCK)])
    query = _write_fasta(tmp_path, "query.fa", [
        ("only_met", "ATG" * 5 + "TAA"),
    ])
    CodonAdaptationIndex(_args(query, ref)).run()
    line = capsys.readouterr().out.strip()
    fields = line.split("\t")
    assert fields[0] == "only_met"
    assert fields[1] == "None"
