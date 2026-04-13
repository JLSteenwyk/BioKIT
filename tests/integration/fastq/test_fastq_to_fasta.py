import pytest

from pathlib import Path
import sys
from unittest.mock import patch

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/DRR284700_1_subset_subset.fq"


@pytest.mark.integration
class TestFastQToFasta(object):
    def test_fastq_to_fasta_stdout(self, capsys):
        testargs = ["biokit", "fastq_to_fasta", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.startswith(">")
        # Should contain FASTA headers only, no '+' quality separator lines
        for line in out.strip().split("\n"):
            assert not line.startswith("+")

    def test_fastq_to_fasta_alias(self, capsys):
        testargs = ["biokit", "fq2fa", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.startswith(">")

    def test_fastq_to_fasta_output_file(self, tmp_path):
        out_path = str(tmp_path / "converted.fa")
        testargs = ["biokit", "fastq_to_fasta", SAMPLE, "-o", out_path]
        with patch.object(sys, "argv", testargs):
            Biokit()
        content = open(out_path).read()
        assert content.startswith(">")
        assert content.count(">") > 0
