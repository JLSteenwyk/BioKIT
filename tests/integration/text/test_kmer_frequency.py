import json

import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/test.nucl.fna"


@pytest.mark.integration
class TestKmerFrequency(object):
    def test_kmer_frequency_basic(self, capsys):
        testargs = ["biokit", "kmer_frequency", SAMPLE, "-k", "3"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert lines[0] == "kmer\tcount"
        assert len(lines) > 1

    def test_kmer_frequency_alias(self, capsys):
        testargs = ["biokit", "kmer_freq", SAMPLE, "-k", "4"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.startswith("kmer\tcount")

    def test_kmer_frequency_canonical(self, capsys):
        testargs = [
            "biokit", "kmer_frequency", SAMPLE, "-k", "3", "--canonical",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        kmers = [line.split("\t")[0] for line in lines[1:]]
        # Canonical should not contain both a kmer and its reverse complement
        # (every reported kmer must be <= its rc lexicographically)
        complement = str.maketrans("ACGT", "TGCA")
        for k in kmers:
            rc = k.translate(complement)[::-1]
            assert k <= rc

    def test_kmer_frequency_verbose(self, capsys):
        testargs = [
            "biokit", "kmer_frequency", SAMPLE, "-k", "3", "-v",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert lines[0] == "id\tkmer\tcount"
        assert len(lines) > 1

    def test_kmer_frequency_json(self, capsys):
        testargs = [
            "biokit", "kmer_frequency", SAMPLE, "-k", "3", "-f", "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        assert len(data) > 0
        for row in data:
            assert set(row.keys()) == {"kmer", "count"}
