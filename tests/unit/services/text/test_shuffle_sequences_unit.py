from argparse import Namespace
from collections import Counter

from biokit.services.text.shuffle_sequences import ShuffleSequences


def test_shuffle_preserves_composition(tmp_path):
    fasta = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    fasta.write_text(">seq1\nAAACCCGGGTTT\n")
    args = Namespace(fasta=str(fasta), output=str(out), seed=42)

    ShuffleSequences(args).run()
    text = out.read_text()
    # Extract sequence (skip header line)
    seq = "".join(line for line in text.strip().split("\n") if not line.startswith(">"))
    assert Counter(seq) == Counter("AAACCCGGGTTT")


def test_shuffle_changes_order_with_seed(tmp_path):
    fasta = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    fasta.write_text(">seq1\nABCDEFGHIJKLMNOP\n")
    args = Namespace(fasta=str(fasta), output=str(out), seed=42)

    ShuffleSequences(args).run()
    text = out.read_text()
    seq = "".join(line for line in text.strip().split("\n") if not line.startswith(">"))
    # With 16 unique chars, shuffling should change the order
    assert seq != "ABCDEFGHIJKLMNOP"


def test_shuffle_reproducible_with_seed(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">seq1\nACGTACGTACGT\n")

    out1 = tmp_path / "out1.fa"
    args1 = Namespace(fasta=str(fasta), output=str(out1), seed=99)
    ShuffleSequences(args1).run()

    out2 = tmp_path / "out2.fa"
    args2 = Namespace(fasta=str(fasta), output=str(out2), seed=99)
    ShuffleSequences(args2).run()

    assert out1.read_text() == out2.read_text()


def test_shuffle_default_output_path(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">seq1\nACGT\n")
    args = Namespace(fasta=str(fasta), output=None, seed=None)

    ShuffleSequences(args).run()
    default_out = tmp_path / "input.fa.shuffled.fa"
    assert default_out.exists()


def test_shuffle_multiple_sequences(tmp_path):
    fasta = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    fasta.write_text(">s1\nAAAA\n>s2\nCCCC\n>s3\nGGGG\n")
    args = Namespace(fasta=str(fasta), output=str(out), seed=42)

    ShuffleSequences(args).run()
    text = out.read_text()
    assert text.count(">") == 3


def test_shuffle_preserves_headers(tmp_path):
    fasta = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    fasta.write_text(">gene1 some description\nACGT\n")
    args = Namespace(fasta=str(fasta), output=str(out), seed=None)

    ShuffleSequences(args).run()
    text = out.read_text()
    assert ">gene1 some description" in text
