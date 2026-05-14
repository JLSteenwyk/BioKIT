import json
import os
from argparse import Namespace

from biokit.services.coding_sequences.effective_number_of_codons import (
    EffectiveNumberOfCodons,
)


def _args(fasta, translation_table=None, verbose=False, plot=False,
          output=None, format=None):
    return Namespace(
        fasta=fasta,
        translation_table=translation_table,
        verbose=verbose,
        plot=plot,
        output=output,
        format=format,
    )


def _write_cds(tmp_path, seqs):
    fa = tmp_path / "cds.fa"
    with open(fa, "w") as fh:
        for name, s in seqs:
            fh.write(f">{name}\n{s}\n")
    return str(fa)


def test_enc_minimum_for_extreme_bias(tmp_path, capsys):
    # Met (1 codon) + Ala (use one of 4 Ala codons exclusively)
    # ENC should be 20: 2 singleton AAs + 18/(F=1)... but only 2 AAs are used.
    # Wright's formula gives the theoretical minimum of 20.
    seq = ("ATG" * 50) + ("GCT" * 50)
    fa = _write_cds(tmp_path, [("biased", seq)])
    EffectiveNumberOfCodons(_args(fa)).run()
    line = capsys.readouterr().out.strip()
    gene_id, enc = line.split("\t")
    assert gene_id == "biased"
    assert float(enc) == 20.0


def test_enc_maximum_for_balanced_six_codon_family(tmp_path, capsys):
    # All 6 Leu codons used equally → F_6 = 1/6 → ENC = 61
    seq = ""
    for codon in ("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"):
        seq += codon * 10
    fa = _write_cds(tmp_path, [("balanced", seq)])
    EffectiveNumberOfCodons(_args(fa)).run()
    line = capsys.readouterr().out.strip()
    _, enc = line.split("\t")
    assert float(enc) == 61.0


def test_enc_skips_sequences_not_divisible_by_3(tmp_path, capsys):
    fa = _write_cds(tmp_path, [
        ("good", "ATGGCTGCT"),  # 9 nt, divisible by 3
        ("bad", "ATGGCTG"),     # 7 nt, not divisible by 3
    ])
    EffectiveNumberOfCodons(_args(fa)).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    assert len(lines) == 1
    assert lines[0].startswith("good")


def test_enc_verbose_includes_gc3(tmp_path, capsys):
    # 100% GC at third position: gene = ATG (Met) + repeated GCG (Ala)
    seq = "ATG" + "GCG" * 30
    fa = _write_cds(tmp_path, [("gc3_high", seq)])
    EffectiveNumberOfCodons(_args(fa, verbose=True)).run()
    fields = capsys.readouterr().out.strip().split("\t")
    assert fields[0] == "gc3_high"
    # ENC parses
    float(fields[1])
    # GC3 should be very high (G at every 3rd position; M's 3rd pos is G too)
    assert float(fields[2]) > 0.9


def test_enc_multiple_genes(tmp_path, capsys):
    fa = _write_cds(tmp_path, [
        ("g1", "ATG" * 30),
        ("g2", "GCT" * 30),
        ("g3", "ATG" * 15 + "GCT" * 15),
    ])
    EffectiveNumberOfCodons(_args(fa)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 3
    ids = [line.split("\t")[0] for line in lines]
    assert ids == ["g1", "g2", "g3"]


def test_enc_json(tmp_path, capsys):
    fa = _write_cds(tmp_path, [
        ("g1", "ATG" * 30),
        ("g2", "GCT" * 30),
    ])
    EffectiveNumberOfCodons(_args(fa, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert len(data) == 2
    for row in data:
        assert set(row.keys()) == {"gene_id", "ENC"}


def test_enc_json_verbose(tmp_path, capsys):
    fa = _write_cds(tmp_path, [("g1", "ATG" * 30)])
    EffectiveNumberOfCodons(_args(fa, verbose=True, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert set(data[0].keys()) == {"gene_id", "ENC", "GC3"}


def test_enc_plot_output(tmp_path):
    fa = _write_cds(tmp_path, [
        ("g1", "ATG" * 30),
        ("g2", "GCT" * 30),
        ("g3", "CTG" * 20 + "TTA" * 10),
    ])
    out_path = str(tmp_path / "enc.png")
    EffectiveNumberOfCodons(_args(fa, plot=True, output=out_path)).run()
    assert os.path.exists(out_path)
    assert os.path.getsize(out_path) > 0


def test_enc_custom_translation_table(tmp_path, capsys):
    # Use bacterial code (11) — should still work
    seq = "ATG" + "GCG" * 20 + "TAA"
    fa = _write_cds(tmp_path, [("bact", seq)])
    EffectiveNumberOfCodons(_args(fa, translation_table="11")).run()
    out = capsys.readouterr().out.strip()
    # One gene reported
    assert out.startswith("bact\t")


def test_enc_expected_curve_at_known_points():
    # Wright's expected ENC: at GC3=0.5, ENC = 2 + 0.5 + 29/(0.25+0.25) = 60.5
    # at GC3=0 or 1, denom=1, ENC = 2 + s + 29 = 31 + s
    assert abs(EffectiveNumberOfCodons._expected_enc(0.5) - 60.5) < 0.01
    assert abs(EffectiveNumberOfCodons._expected_enc(0.0) - 31.0) < 0.01
    assert abs(EffectiveNumberOfCodons._expected_enc(1.0) - 32.0) < 0.01
