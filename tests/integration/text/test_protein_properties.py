import json

import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/protein.fa"


@pytest.mark.integration
class TestProteinProperties(object):
    def test_protein_properties_basic(self, capsys):
        testargs = ["biokit", "protein_properties", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        header = lines[0].split("\t")
        assert header[0] == "id"
        assert "molecular_weight" in header
        assert "isoelectric_point" in header
        assert "gravy" in header
        assert "aromaticity" in header
        assert "instability_index" in header
        # protein.fa has 3 sequences
        assert len(lines) == 4

    def test_protein_properties_alias(self, capsys):
        testargs = ["biokit", "prot_prop", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.startswith("id\t")

    def test_protein_properties_json(self, capsys):
        testargs = ["biokit", "protein_properties", SAMPLE, "-f", "json"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        assert len(data) == 3
        for row in data:
            assert {
                "id", "length", "molecular_weight", "isoelectric_point",
                "gravy", "aromaticity", "instability_index",
            } == set(row.keys())
