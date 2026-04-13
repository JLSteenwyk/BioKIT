import os

import pytest

from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/test_genbank.gb"


@pytest.mark.integration
class TestGenbankToFasta(object):
    def test_genbank_to_fasta_whole_record(self, capsys):
        testargs = ["biokit", "genbank_to_fasta", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert ">TEST001.1" in out
        assert out.count(">") == 1

    def test_genbank_to_fasta_alias(self, capsys):
        testargs = ["biokit", "gb2fa", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert ">TEST001.1" in out

    def test_genbank_to_fasta_cds_filter(self, capsys):
        testargs = ["biokit", "genbank_to_fasta", SAMPLE, "-t", "CDS"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.count(">") == 2
        assert ">TEST001_01" in out
        assert ">TEST001_02" in out

    def test_genbank_to_fasta_translate(self, capsys):
        testargs = [
            "biokit",
            "genbank_to_fasta",
            SAMPLE,
            "-t",
            "CDS",
            "--translate",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert "MKIL*" in out
        assert "MARN*" in out

    def test_genbank_to_fasta_output_file(self, tmp_path):
        out_path = str(tmp_path / "out.fa")
        testargs = [
            "biokit",
            "genbank_to_fasta",
            SAMPLE,
            "-t",
            "rRNA",
            "-o",
            out_path,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert os.path.exists(out_path)
        content = open(out_path).read()
        assert ">TEST001_rrna_01" in content
