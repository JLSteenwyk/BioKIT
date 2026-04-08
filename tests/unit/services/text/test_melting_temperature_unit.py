from argparse import Namespace

from biokit.services.text.melting_temperature import MeltingTemperature


def test_tm_default_conditions(tmp_path, capsys):
    fasta = tmp_path / "primers.fa"
    fasta.write_text(">p1\nATGCATGCATGC\n>p2\nGCGCGCGCGCGC\n")
    args = Namespace(fasta=str(fasta), na=50.0, oligo_conc=250.0, format=None)

    MeltingTemperature(args).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 2
    # GC-rich sequence should have higher Tm
    tm_p1 = float(lines[0].split("\t")[1])
    tm_p2 = float(lines[1].split("\t")[1])
    assert tm_p2 > tm_p1


def test_tm_higher_na_raises_tm(tmp_path, capsys):
    fasta = tmp_path / "primer.fa"
    fasta.write_text(">p1\nATGCATGCATGC\n")

    args_low = Namespace(fasta=str(fasta), na=10.0, oligo_conc=250.0, format=None)
    MeltingTemperature(args_low).run()
    tm_low = float(capsys.readouterr().out.strip().split("\t")[1])

    args_high = Namespace(fasta=str(fasta), na=200.0, oligo_conc=250.0, format=None)
    MeltingTemperature(args_high).run()
    tm_high = float(capsys.readouterr().out.strip().split("\t")[1])

    assert tm_high > tm_low


def test_tm_defaults(tmp_path):
    fasta = tmp_path / "primer.fa"
    fasta.write_text(">p1\nATGC\n")
    args = Namespace(fasta=str(fasta), na=None, oligo_conc=None, format=None)

    svc = MeltingTemperature(args)
    assert svc.na_conc == 50.0
    assert svc.oligo_conc == 250.0


def test_tm_json_format(tmp_path, capsys):
    fasta = tmp_path / "primer.fa"
    fasta.write_text(">p1\nATGCATGC\n")
    args = Namespace(fasta=str(fasta), na=50.0, oligo_conc=250.0, format="json")

    MeltingTemperature(args).run()
    import json
    out = json.loads(capsys.readouterr().out)
    assert len(out) == 1
    assert out[0]["id"] == "p1"
    assert "Tm" in out[0]


def test_tm_single_sequence(tmp_path, capsys):
    fasta = tmp_path / "primer.fa"
    fasta.write_text(">oligo\nGATCGATCGATCGATC\n")
    args = Namespace(fasta=str(fasta), na=50.0, oligo_conc=250.0, format=None)

    MeltingTemperature(args).run()
    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 1
    assert lines[0].startswith("oligo\t")
    tm = float(lines[0].split("\t")[1])
    assert 30 < tm < 80
