import json

import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/test.nucl.fna"


@pytest.mark.integration
class TestFindOrfs(object):
    def test_find_orfs_basic(self, capsys):
        testargs = [
            "biokit", "find_orfs", SAMPLE, "-m", "30",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert lines[0] == (
            "id\torf_id\tframe\tstart\tstop\tlength_nt\tlength_aa"
        )

    def test_find_orfs_alias(self, capsys):
        testargs = ["biokit", "orfs", SAMPLE, "-m", "30"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.startswith(
            "id\torf_id\tframe\tstart\tstop\tlength_nt\tlength_aa"
        )

    def test_find_orfs_extract(self, capsys):
        testargs = [
            "biokit", "find_orfs", SAMPLE, "-m", "30", "--extract",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        # When --extract is used, output should be FASTA (every other line
        # is a header beginning with '>')
        lines = [line for line in out.strip().split("\n") if line]
        if lines:
            headers = [ln for ln in lines if ln.startswith(">")]
            sequences = [ln for ln in lines if not ln.startswith(">")]
            assert len(headers) == len(sequences)
            # Each nucleotide ORF should begin with ATG
            for seq in sequences:
                assert seq.startswith("ATG")

    def test_find_orfs_extract_protein(self, capsys):
        testargs = [
            "biokit", "find_orfs", SAMPLE, "-m", "30",
            "--extract", "--protein",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = [line for line in out.strip().split("\n") if line]
        sequences = [ln for ln in lines if not ln.startswith(">")]
        for seq in sequences:
            assert seq.startswith("M")
            assert "*" not in seq

    def test_find_orfs_json(self, capsys):
        testargs = [
            "biokit", "find_orfs", SAMPLE, "-m", "30", "-f", "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        for row in data:
            assert set(row.keys()) == {
                "id", "orf_id", "frame", "start", "stop",
                "length_nt", "length_aa",
            }
