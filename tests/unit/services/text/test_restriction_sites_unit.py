from argparse import Namespace

from biokit.services.text.restriction_sites import RestrictionSites


def test_single_enzyme_with_sites(tmp_path, capsys):
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">s1\nGAATTCAAAAAAGAATTC\n")
    args = Namespace(fasta=str(fasta), enzymes=["EcoRI"])

    RestrictionSites(args).run()
    out = capsys.readouterr().out.strip()
    parts = out.split("\t")
    assert parts[0] == "s1"
    assert parts[1] == "EcoRI"
    assert parts[2] == "2"


def test_enzyme_no_sites(tmp_path, capsys):
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">s1\nAAAAAAAA\n")
    args = Namespace(fasta=str(fasta), enzymes=["EcoRI"])

    RestrictionSites(args).run()
    out = capsys.readouterr().out.strip()
    assert "no_sites" in out


def test_multiple_enzymes(tmp_path, capsys):
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">s1\nGAATTCAAGGATCC\n")
    args = Namespace(fasta=str(fasta), enzymes=["EcoRI,BamHI"])

    RestrictionSites(args).run()
    out = capsys.readouterr().out.strip()
    lines = out.split("\n")
    assert len(lines) == 2
    enzymes_found = {line.split("\t")[1] for line in lines}
    assert enzymes_found == {"EcoRI", "BamHI"}


def test_fragment_sizes(tmp_path, capsys):
    # GAATTC at position 1, cuts after pos 1
    # Sequence: GAATTCGAATTC (len=12), EcoRI sites at 2 and 8
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">s1\nGAATTCGAATTC\n")
    args = Namespace(fasta=str(fasta), enzymes=["EcoRI"])

    RestrictionSites(args).run()
    out = capsys.readouterr().out.strip()
    parts = out.split("\t")
    # Should have fragment sizes in last column
    fragments = [int(f) for f in parts[4].split(",")]
    assert sum(fragments) == 12


def test_multiple_sequences(tmp_path, capsys):
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">s1\nGAATTC\n>s2\nAAAA\n")
    args = Namespace(fasta=str(fasta), enzymes=["EcoRI"])

    RestrictionSites(args).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 2
    assert lines[0].startswith("s1")
    assert lines[1].startswith("s2")


def test_comma_separated_enzymes_arg(tmp_path, capsys):
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">s1\nGAATTCGGATCC\n")
    args = Namespace(fasta=str(fasta), enzymes=["EcoRI", "BamHI"])

    RestrictionSites(args).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 2
