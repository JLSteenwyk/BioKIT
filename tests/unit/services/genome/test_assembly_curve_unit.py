import json
import os
from argparse import Namespace

from biokit.services.genome.assembly_curve import AssemblyCurve


def test_assembly_curve_tsv(tmp_path, capsys):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nATCGATCG\n>c2\nATCG\n>c3\nATCGATCGATCG\n")
    args = Namespace(fasta=str(fasta), format=None, plot=False, output=None)
    AssemblyCurve(args).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    assert lines[0] == "rank\tlength\tcumulative_length"
    assert lines[1] == "1\t12\t12"
    assert lines[2] == "2\t8\t20"
    assert lines[3] == "3\t4\t24"


def test_assembly_curve_json(tmp_path, capsys):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nATCGATCG\n>c2\nATCG\n>c3\nATCGATCGATCG\n")
    args = Namespace(fasta=str(fasta), format="json", plot=False, output=None)
    AssemblyCurve(args).run()
    out = capsys.readouterr().out
    data = json.loads(out)
    assert len(data) == 3
    assert data[0]["rank"] == 1
    assert data[0]["length"] == 12
    assert data[0]["cumulative_length"] == 12
    assert data[2]["cumulative_length"] == 24


def test_assembly_curve_yaml(tmp_path, capsys):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nATCGATCG\n>c2\nATCG\n>c3\nATCGATCGATCG\n")
    args = Namespace(fasta=str(fasta), format="yaml", plot=False, output=None)
    AssemblyCurve(args).run()
    out = capsys.readouterr().out
    assert "rank: 1" in out
    assert "cumulative_length: 24" in out


def test_assembly_curve_plot(tmp_path, capsys):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nATCGATCG\n>c2\nATCG\n>c3\nATCGATCGATCG\n")
    plot_path = str(tmp_path / "curve.png")
    args = Namespace(fasta=str(fasta), format=None, plot=True, output=plot_path)
    AssemblyCurve(args).run()
    assert os.path.exists(plot_path)
    assert os.path.getsize(plot_path) > 0


def test_assembly_curve_single_contig(tmp_path, capsys):
    fasta = tmp_path / "single.fa"
    fasta.write_text(">c1\nATCGATCG\n")
    args = Namespace(fasta=str(fasta), format=None, plot=False, output=None)
    AssemblyCurve(args).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    assert len(lines) == 2
    assert lines[1] == "1\t8\t8"


def test_assembly_curve_sorted_descending(tmp_path, capsys):
    fasta = tmp_path / "contigs.fa"
    # Write in non-sorted order; output should still be descending
    fasta.write_text(">small\nAT\n>big\nATCGATCGATCG\n>med\nATCGAT\n")
    args = Namespace(fasta=str(fasta), format=None, plot=False, output=None)
    AssemblyCurve(args).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    lengths = [int(line.split("\t")[1]) for line in lines[1:]]
    assert lengths == sorted(lengths, reverse=True)
