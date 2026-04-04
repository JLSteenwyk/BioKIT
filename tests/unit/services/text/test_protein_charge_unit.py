from argparse import Namespace

from biokit.services.text.protein_charge import ProteinCharge


def test_charge_at_default_pH(tmp_path, capsys):
    fasta = tmp_path / "prot.fa"
    fasta.write_text(">seq1\nKKKK\n>seq2\nDDDD\n")
    args = Namespace(fasta=str(fasta), pH=7.0, format=None)

    ProteinCharge(args).run()
    out = capsys.readouterr().out.strip().split("\n")
    assert len(out) == 2
    # KKKK is highly positive at pH 7
    assert out[0].startswith("seq1\t")
    charge1 = float(out[0].split("\t")[1])
    assert charge1 > 3.0
    # DDDD is negative at pH 7
    assert out[1].startswith("seq2\t")
    charge2 = float(out[1].split("\t")[1])
    assert charge2 < -3.0


def test_charge_at_custom_pH(tmp_path, capsys):
    fasta = tmp_path / "prot.fa"
    fasta.write_text(">seq1\nKKKK\n")
    args = Namespace(fasta=str(fasta), pH=2.0, format=None)

    ProteinCharge(args).run()
    out = capsys.readouterr().out.strip()
    charge = float(out.split("\t")[1])
    # At very low pH, lysines are fully protonated -> higher positive charge
    assert charge > 3.5


def test_charge_default_pH_is_7(tmp_path, capsys):
    fasta = tmp_path / "prot.fa"
    fasta.write_text(">seq1\nACDE\n")
    args = Namespace(fasta=str(fasta), pH=None, format=None)

    svc = ProteinCharge(args)
    assert svc.pH == 7.0


def test_charge_single_sequence(tmp_path, capsys):
    fasta = tmp_path / "prot.fa"
    fasta.write_text(">only\nMKWVTFISLLFLFSSAYS\n")
    args = Namespace(fasta=str(fasta), pH=7.0, format=None)

    ProteinCharge(args).run()
    out = capsys.readouterr().out.strip().split("\n")
    assert len(out) == 1
    assert out[0].startswith("only\t")


def test_charge_json_format(tmp_path, capsys):
    fasta = tmp_path / "prot.fa"
    fasta.write_text(">seq1\nKKKK\n")
    args = Namespace(fasta=str(fasta), pH=7.0, format="json")

    ProteinCharge(args).run()
    import json
    out = json.loads(capsys.readouterr().out)
    assert len(out) == 1
    assert out[0]["id"] == "seq1"
    assert "charge" in out[0]
    assert out[0]["pH"] == 7.0
