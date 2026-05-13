import json
from argparse import Namespace

from biokit.services.text.protein_properties import (
    PROPERTY_FIELDS,
    ProteinProperties,
)


def _args(fasta, format=None):
    return Namespace(fasta=fasta, format=format)


def _write_fasta(tmp_path, seqs):
    fa = tmp_path / "prot.fa"
    with open(fa, "w") as fh:
        for name, seq in seqs:
            fh.write(f">{name}\n{seq}\n")
    return str(fa)


def test_protein_properties_tsv_header_and_columns(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("a", "MKTAYIAKQRQISFVKSHFSRQ")])
    ProteinProperties(_args(fa)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0] == "\t".join(PROPERTY_FIELDS)
    assert len(lines) == 2
    fields = lines[1].split("\t")
    assert fields[0] == "a"
    # length is amino acids
    assert fields[1] == "22"
    # MW, pI, gravy, aromaticity, instability_index are all numeric
    for value in fields[2:]:
        float(value)  # should parse


def test_protein_properties_multiple_sequences(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [
        ("a", "MKTAYIAKQR"),
        ("b", "GIVEQCCTSI"),
        ("c", "WWWFYY"),
    ])
    ProteinProperties(_args(fa)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 4  # header + 3 rows
    ids = [line.split("\t")[0] for line in lines[1:]]
    assert ids == ["a", "b", "c"]


def test_protein_properties_aromaticity(tmp_path, capsys):
    # All aromatic (F, W, Y) → aromaticity = 1.0
    fa = _write_fasta(tmp_path, [("aromatic", "WWWWYYYYFFFF")])
    ProteinProperties(_args(fa)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    fields = lines[1].split("\t")
    aromaticity = float(fields[PROPERTY_FIELDS.index("aromaticity")])
    assert aromaticity == 1.0

    # No aromatic → aromaticity = 0
    fa2 = _write_fasta(tmp_path, [("non_aromatic", "AAAAAAAAAAAA")])
    ProteinProperties(_args(fa2)).run()
    lines2 = capsys.readouterr().out.strip().split("\n")
    fields2 = lines2[1].split("\t")
    aromaticity2 = float(fields2[PROPERTY_FIELDS.index("aromaticity")])
    assert aromaticity2 == 0.0


def test_protein_properties_strips_stop_and_gaps(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("with_stop", "MKTAY*"), ("with_gap", "MK-TAY")])
    ProteinProperties(_args(fa)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    # Both sequences should have length 5 after stripping
    a_fields = lines[1].split("\t")
    b_fields = lines[2].split("\t")
    assert a_fields[1] == "5"
    assert b_fields[1] == "5"


def test_protein_properties_case_insensitive(tmp_path, capsys):
    fa1 = _write_fasta(tmp_path, [("a", "mktayiak")])
    ProteinProperties(_args(fa1)).run()
    lower_out = capsys.readouterr().out

    fa2 = _write_fasta(tmp_path, [("a", "MKTAYIAK")])
    ProteinProperties(_args(fa2)).run()
    upper_out = capsys.readouterr().out

    assert lower_out == upper_out


def test_protein_properties_json(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [
        ("a", "MKTAYIAK"),
        ("b", "GIVEQCCT"),
    ])
    ProteinProperties(_args(fa, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert len(data) == 2
    for row in data:
        assert set(row.keys()) == set(PROPERTY_FIELDS)


def test_protein_properties_yaml(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("a", "MKTAYIAK")])
    ProteinProperties(_args(fa, format="yaml")).run()
    out = capsys.readouterr().out
    assert "id: a" in out
    assert "molecular_weight:" in out
    assert "isoelectric_point:" in out
    assert "gravy:" in out


def test_protein_properties_empty_sequence(tmp_path, capsys):
    fa = _write_fasta(tmp_path, [("empty", "")])
    ProteinProperties(_args(fa)).run()
    lines = capsys.readouterr().out.strip().split("\n")
    fields = lines[1].split("\t")
    assert fields[0] == "empty"
    assert fields[1] == "0"
    # All properties should be "None"
    for value in fields[2:]:
        assert value == "None"


def test_protein_properties_known_values(tmp_path, capsys):
    # A small peptide. ProtParam gives reproducible numbers.
    fa = _write_fasta(tmp_path, [("p", "MFFFY")])
    ProteinProperties(_args(fa)).run()
    fields = capsys.readouterr().out.strip().split("\n")[1].split("\t")
    # length is 5
    assert fields[1] == "5"
    mw = float(fields[PROPERTY_FIELDS.index("molecular_weight")])
    # MFFFY: M (149.21) + F (165.19)*3 + Y (181.19) − 4 waters (4*18.015)
    # = 149.21 + 495.57 + 181.19 − 72.06 = 753.91
    assert abs(mw - 753.91) < 0.5
    aromaticity = float(fields[PROPERTY_FIELDS.index("aromaticity")])
    # 4 of 5 are aromatic (F, F, F, Y)
    assert aromaticity == 0.8
