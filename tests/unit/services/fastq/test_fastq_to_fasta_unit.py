from argparse import Namespace
from io import StringIO
from unittest.mock import patch

from biokit.services.fastq.fastq_to_fasta import FastQToFasta


SAMPLE_FASTQ = (
    "@read1\n"
    "ATCGATCG\n"
    "+\n"
    "IIIIIIII\n"
    "@read2 description\n"
    "GGGGCCCC\n"
    "+\n"
    "HHHHHHHH\n"
    "@read3\n"
    "NNNNNNNN\n"
    "+\n"
    "!!!!!!!!\n"
)


def _args(fastq, output=None):
    return Namespace(fastq=fastq, output=output)


def test_fastq_to_fasta_stdout(tmp_path, capsys):
    fq = tmp_path / "reads.fq"
    fq.write_text(SAMPLE_FASTQ)
    FastQToFasta(_args(str(fq))).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    assert lines[0] == ">read1"
    assert lines[1] == "ATCGATCG"
    assert lines[2] == ">read2 description"
    assert lines[3] == "GGGGCCCC"
    assert lines[4] == ">read3"
    assert lines[5] == "NNNNNNNN"


def test_fastq_to_fasta_output_file(tmp_path):
    fq = tmp_path / "reads.fq"
    fq.write_text(SAMPLE_FASTQ)
    out_path = tmp_path / "out.fa"
    FastQToFasta(_args(str(fq), output=str(out_path))).run()
    content = out_path.read_text()
    assert content.count(">") == 3
    assert ">read1\nATCGATCG" in content
    assert ">read2 description\nGGGGCCCC" in content


def test_fastq_to_fasta_record_count(tmp_path, capsys):
    fq = tmp_path / "reads.fq"
    fq.write_text(SAMPLE_FASTQ)
    FastQToFasta(_args(str(fq))).run()
    out = capsys.readouterr().out
    assert out.count(">") == 3
    # No quality score lines should leak through
    assert "+\n" not in out
    assert "IIII" not in out


def test_fastq_to_fasta_stdin(capsys):
    with patch("sys.stdin", StringIO(SAMPLE_FASTQ)):
        FastQToFasta(_args("-")).run()
    out = capsys.readouterr().out
    assert out.count(">") == 3
    assert ">read1" in out


def test_fastq_to_fasta_empty_file(tmp_path, capsys):
    fq = tmp_path / "empty.fq"
    fq.write_text("")
    FastQToFasta(_args(str(fq))).run()
    out = capsys.readouterr().out
    assert out == ""
