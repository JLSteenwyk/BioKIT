from argparse import Namespace

from biokit.services.text.fasta_deduplication import FastaDeduplication


def test_dedup_removes_exact_duplicates(tmp_path):
    fasta = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    fasta.write_text(">a\nACGT\n>b\nTTTT\n>c\nACGT\n")
    args = Namespace(fasta=str(fasta), output=str(out))

    FastaDeduplication(args).run()
    text = out.read_text()
    assert text.count(">") == 2
    assert ">a" in text
    assert ">b" in text
    assert ">c" not in text


def test_dedup_case_insensitive(tmp_path):
    fasta = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    fasta.write_text(">a\nACGT\n>b\nacgt\n")
    args = Namespace(fasta=str(fasta), output=str(out))

    FastaDeduplication(args).run()
    text = out.read_text()
    assert text.count(">") == 1
    assert ">a" in text


def test_dedup_no_duplicates(tmp_path):
    fasta = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    fasta.write_text(">a\nACGT\n>b\nTTTT\n>c\nGGGG\n")
    args = Namespace(fasta=str(fasta), output=str(out))

    FastaDeduplication(args).run()
    text = out.read_text()
    assert text.count(">") == 3


def test_dedup_default_output_path(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">a\nACGT\n>b\nACGT\n")
    args = Namespace(fasta=str(fasta), output=None)

    FastaDeduplication(args).run()
    default_out = tmp_path / "input.fa.dedup.fa"
    assert default_out.exists()
    text = default_out.read_text()
    assert text.count(">") == 1


def test_dedup_single_sequence(tmp_path):
    fasta = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    fasta.write_text(">a\nACGT\n")
    args = Namespace(fasta=str(fasta), output=str(out))

    FastaDeduplication(args).run()
    text = out.read_text()
    assert text.count(">") == 1
    assert ">a" in text
