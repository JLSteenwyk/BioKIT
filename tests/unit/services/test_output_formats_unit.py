import json
from argparse import Namespace

import pytest

from biokit.services.alignment.alignment_summary import AlignmentSummary
from biokit.services.alignment.constant_sites import ConstantSites
from biokit.services.coding_sequences.relative_synonymous_codon_usage import (
    RelativeSynonymousCodonUsage,
)
from biokit.services.fastq.fastq_read_lengths import FastQReadLengths
from biokit.services.fastq.trim_pe_fastq import TrimPEFastQ
from biokit.services.fastq.trim_se_fastq import TrimSEFastQ
from biokit.services.genome.genome_assembly_metrics import GenomeAssemblyMetrics
from biokit.services.genome.gc_content import GCContent
from biokit.services.text.character_frequency import CharacterFrequency
from biokit.services.text.sequence_length import SequenceLength


def test_character_frequency_json_output(tmp_path, capsys):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">a\nACGT\n")
    args = Namespace(fasta=str(fasta), format="json")

    CharacterFrequency(args).run()
    out = capsys.readouterr().out.strip()
    data = json.loads(out)
    assert data == [
        {"character": "A", "frequency": 0.25},
        {"character": "C", "frequency": 0.25},
        {"character": "G", "frequency": 0.25},
        {"character": "T", "frequency": 0.25},
    ]


def test_sequence_length_yaml_output_is_sorted(tmp_path, capsys):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">z\nAA\n>a\nAAAA\n")
    args = Namespace(fasta=str(fasta), format="yaml")

    SequenceLength(args).run()
    out = capsys.readouterr().out.strip().splitlines()
    assert out[0] == "- id: a"
    assert "  length: 4" in out
    assert "- id: z" in out


def test_rscu_json_output_with_custom_table(tmp_path, capsys):
    fasta = tmp_path / "cds.fa"
    table = tmp_path / "tt.txt"
    fasta.write_text(">g1\nATGTTTTTC\n")
    table.write_text("AUG M\nUUU F\nUUC F\n")
    args = Namespace(fasta=str(fasta), translation_table=str(table), format="json")

    RelativeSynonymousCodonUsage(args).run()
    data = json.loads(capsys.readouterr().out.strip())
    assert data == [
        {"codon": "AUG", "rscu": 1.0},
        {"codon": "UUC", "rscu": 1.0},
        {"codon": "UUU", "rscu": 1.0},
    ]


def test_genome_assembly_metrics_json_output(tmp_path, capsys):
    fasta = tmp_path / "g.fa"
    fasta.write_text(">a\nATGC\n")
    args = Namespace(fasta=str(fasta), threshold=None, format="json")

    GenomeAssemblyMetrics(args).run()
    data = json.loads(capsys.readouterr().out.strip())
    assert data["assembly_size"] == 4
    assert data["gc_content"] == 0.5
    assert data["number_of_scaffolds"] == 1


def test_invalid_output_format_raises_value_error(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">a\nAC\n")
    args = Namespace(fasta=str(fasta), format="xml")

    with pytest.raises(ValueError, match="Invalid output format"):
        SequenceLength(args).run()


def test_alignment_summary_json_output(tmp_path, capsys):
    fasta = tmp_path / "aln.fa"
    fasta.write_text(">a\nA-\n>b\nAA\n")
    args = Namespace(fasta=str(fasta), format="json")

    AlignmentSummary(args).run()
    data = json.loads(capsys.readouterr().out.strip())
    assert data["number_of_taxa"] == 2
    assert data["alignment_length"] == 2
    assert "character_frequencies" in data


def test_constant_sites_verbose_yaml_output(tmp_path, capsys):
    fasta = tmp_path / "aln.fa"
    fasta.write_text(">a\nAA\n>b\nAT\n>c\nAA\n")
    args = Namespace(fasta=str(fasta), verbose=True, format="yaml")

    ConstantSites(args).run()
    out = capsys.readouterr().out
    assert "- site_index: 0" in out
    assert "classification: constant_site" in out


def test_gc_content_verbose_json_output(tmp_path, capsys):
    fasta = tmp_path / "g.fa"
    fasta.write_text(">z\nGG\n>a\nAT\n")
    args = Namespace(fasta=str(fasta), verbose=True, format="json")

    GCContent(args).run()
    data = json.loads(capsys.readouterr().out.strip())
    assert data == [
        {"entry": "a", "gc_content": 0.0},
        {"entry": "z", "gc_content": 1.0},
    ]


def test_fastq_read_lengths_json_output(tmp_path, capsys):
    fastq = tmp_path / "reads.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n@r2\nACG\n+\nIII\n")
    args = Namespace(fastq=str(fastq), verbose=False, format="json")

    FastQReadLengths(args).run()
    data = json.loads(capsys.readouterr().out.strip())
    assert data["mean"] == 3.5


def test_trim_se_fastq_json_summary(tmp_path, capsys):
    fastq = tmp_path / "reads.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGA\n+\nIIII\n")
    out = tmp_path / "trimmed.fq"
    args = Namespace(
        fastq=str(fastq),
        minimum=None,
        length=None,
        output_file=str(out),
        format="json",
    )

    TrimSEFastQ(args).run()
    data = json.loads(capsys.readouterr().out.strip())
    assert data == {"reads_processed": 2, "reads_kept": 2, "reads_removed": 0}


def test_trim_pe_fastq_json_summary(tmp_path, capsys):
    fastq1 = tmp_path / "r1.fq"
    fastq2 = tmp_path / "r2.fq"
    fastq1.write_text("@r1\nACGT\n+\nIIII\n")
    fastq2.write_text("@r1\nTGCA\n+\nIIII\n")
    args = Namespace(
        fastq1=str(fastq1),
        fastq2=str(fastq2),
        minimum=None,
        length=None,
        format="json",
    )

    TrimPEFastQ(args).run()
    data = json.loads(capsys.readouterr().out.strip())
    assert data["reads_processed"] == 2
    assert data["pairs_kept"] == 1
