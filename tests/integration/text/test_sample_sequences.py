import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/test.nucl.fna"


@pytest.mark.integration
class TestSampleSequences(object):
    def test_sample_sequences_number(self, capsys):
        testargs = [
            "biokit", "sample_sequences", SAMPLE, "-n", "2", "-s", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.count(">") == 2

    def test_sample_sequences_alias(self, capsys):
        testargs = [
            "biokit", "sample_seqs", SAMPLE, "-n", "1", "-s", "0",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.count(">") == 1

    def test_sample_sequences_seed_reproducible(self, capsys):
        testargs = [
            "biokit", "sample_sequences", SAMPLE, "-n", "2", "-s", "7",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        first = capsys.readouterr().out

        with patch.object(sys, "argv", testargs):
            Biokit()
        second = capsys.readouterr().out

        assert first == second

    def test_sample_sequences_output_file(self, tmp_path):
        out_path = str(tmp_path / "out.fa")
        testargs = [
            "biokit", "sample_sequences", SAMPLE,
            "-n", "2", "-s", "0", "-o", out_path,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        content = open(out_path).read()
        assert content.count(">") == 2
