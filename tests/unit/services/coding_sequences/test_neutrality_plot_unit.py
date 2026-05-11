import json
import os
from argparse import Namespace

from biokit.services.coding_sequences.neutrality_plot import NeutralityPlot


def _args(fasta, verbose=False, plot=False, output=None, format=None):
    return Namespace(
        fasta=fasta,
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


def test_neutrality_plot_summary_tsv(tmp_path, capsys):
    # 5 CDS with varying GC content
    seqs = [
        ("g1", "ATGAAAAAA" * 10),
        ("g2", "ATGCCCCCC" * 10),
        ("g3", "ATGGGGGGG" * 10),
        ("g4", "ATGTCTGCG" * 10),
        ("g5", "ATGCGCGCG" * 10),
    ]
    fa = _write_cds(tmp_path, seqs)
    NeutralityPlot(_args(fa)).run()
    out = capsys.readouterr().out.strip().split("\n")
    # tsv summary has 4 lines: slope, intercept, r_squared, n
    assert len(out) == 4
    labels = [line.split("\t")[1] for line in out]
    assert labels == ["slope", "intercept", "r_squared", "n"]


def test_neutrality_plot_verbose_tsv(tmp_path, capsys):
    seqs = [
        ("g1", "ATGAAAAAA" * 5),
        ("g2", "ATGCCCCCC" * 5),
        ("g3", "ATGGGGGGG" * 5),
    ]
    fa = _write_cds(tmp_path, seqs)
    NeutralityPlot(_args(fa, verbose=True)).run()
    out = capsys.readouterr().out
    assert "id\tGC12\tGC3" in out
    assert "g1\t" in out
    assert "g2\t" in out
    assert "g3\t" in out
    # Summary lines also present
    assert "slope" in out
    assert "r_squared" in out


def test_neutrality_plot_known_values(tmp_path, capsys):
    # A single gene of pure A's at all three codon positions
    # All GC = 0, regression undefined (gc3 std == 0) → slope=0
    seqs = [
        ("a", "AAAAAAAAA" * 5),
        ("b", "AAAAAAAAA" * 5),
    ]
    fa = _write_cds(tmp_path, seqs)
    NeutralityPlot(_args(fa)).run()
    out = capsys.readouterr().out
    # slope should be 0 since gc3 is constant
    slope = float(out.splitlines()[0].split("\t")[0])
    assert slope == 0.0


def test_neutrality_plot_json(tmp_path, capsys):
    seqs = [
        ("g1", "ATGAAAAAA" * 10),
        ("g2", "ATGCCCCCC" * 10),
        ("g3", "ATGGGGGGG" * 10),
        ("g4", "ATGTCTGCG" * 10),
    ]
    fa = _write_cds(tmp_path, seqs)
    NeutralityPlot(_args(fa, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert "regression" in data
    assert set(data["regression"].keys()) == {"slope", "intercept", "r_squared", "n"}
    assert data["regression"]["n"] == 4
    assert "per_gene" not in data


def test_neutrality_plot_json_verbose(tmp_path, capsys):
    seqs = [
        ("g1", "ATGAAAAAA" * 5),
        ("g2", "ATGCCCCCC" * 5),
    ]
    fa = _write_cds(tmp_path, seqs)
    NeutralityPlot(_args(fa, verbose=True, format="json")).run()
    data = json.loads(capsys.readouterr().out)
    assert "regression" in data
    assert "per_gene" in data
    assert len(data["per_gene"]) == 2
    for row in data["per_gene"]:
        assert set(row.keys()) == {"id", "GC12", "GC3"}


def test_neutrality_plot_plot_output(tmp_path):
    seqs = [
        ("g1", "ATGAAAAAA" * 10),
        ("g2", "ATGCCCCCC" * 10),
        ("g3", "ATGGGGGGG" * 10),
        ("g4", "ATGTCTGCG" * 10),
    ]
    fa = _write_cds(tmp_path, seqs)
    out_path = str(tmp_path / "neut.png")
    NeutralityPlot(_args(fa, plot=True, output=out_path)).run()
    assert os.path.exists(out_path)
    assert os.path.getsize(out_path) > 0


def test_neutrality_plot_single_gene(tmp_path, capsys):
    seqs = [("only", "ATGCGCGCG" * 10)]
    fa = _write_cds(tmp_path, seqs)
    NeutralityPlot(_args(fa)).run()
    out = capsys.readouterr().out
    # With n=1, regression yields None for slope/intercept/r_squared
    lines = out.strip().split("\n")
    n_line = [ln for ln in lines if ln.endswith("\tn")][0]
    assert n_line.startswith("1\t")
    slope_line = [ln for ln in lines if ln.endswith("\tslope")][0]
    assert slope_line.startswith("None\t")
