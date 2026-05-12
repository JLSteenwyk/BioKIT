import json

import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/test.nucl.fna"


@pytest.mark.integration
class TestDinucleotideOdds(object):
    def test_dinucleotide_odds_basic(self, capsys):
        testargs = ["biokit", "dinucleotide_odds", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert lines[0] == "dinucleotide\tO_E"
        assert len(lines) == 17

    def test_dinucleotide_odds_alias(self, capsys):
        testargs = ["biokit", "dno", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.startswith("dinucleotide\tO_E")

    def test_dinucleotide_odds_verbose(self, capsys):
        testargs = ["biokit", "dinucleotide_odds", SAMPLE, "-v"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        header = lines[0].split("\t")
        assert header[0] == "id"
        # All 16 dinucleotides as columns
        assert len(header) == 17

    def test_dinucleotide_odds_json(self, capsys):
        testargs = ["biokit", "dinucleotide_odds", SAMPLE, "-f", "json"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        assert len(data) == 16
        for row in data:
            assert set(row.keys()) == {"dinucleotide", "O_E"}
