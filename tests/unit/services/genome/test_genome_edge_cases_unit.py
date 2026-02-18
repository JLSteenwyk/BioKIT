from argparse import Namespace

from biokit.services.genome.gc_content import GCContent
from biokit.services.genome.genome_assembly_metrics import GenomeAssemblyMetrics


def test_gc_content_handles_gap_only_sequences_non_verbose(tmp_path, capsys):
    fasta = tmp_path / "gap_only.fa"
    fasta.write_text(">a\n---???---\n>b\n??--\n")
    args = Namespace(fasta=str(fasta), verbose=False)

    GCContent(args).run()
    assert capsys.readouterr().out.strip() == "0"


def test_gc_content_handles_gap_only_sequences_verbose(tmp_path, capsys):
    fasta = tmp_path / "gap_only.fa"
    fasta.write_text(">a\n---???---\n")
    args = Namespace(fasta=str(fasta), verbose=True)

    GCContent(args).run()
    assert capsys.readouterr().out.strip() == "a\t0"


def test_genome_assembly_metrics_handles_empty_fasta(tmp_path, capsys):
    fasta = tmp_path / "empty.fa"
    fasta.write_text("")
    args = Namespace(fasta=str(fasta), threshold=None)

    GenomeAssemblyMetrics(args).run()
    out = capsys.readouterr().out
    assert "0\tAssembly size" in out
    assert "0.0\tGC content" in out
