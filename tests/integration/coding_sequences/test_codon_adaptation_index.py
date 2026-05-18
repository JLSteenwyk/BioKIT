import json

import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = (
    f"{here.parent.parent.parent}/sample_files/"
    "GCF_000146045.2_R64_cds_from_genomic.small.fna"
)


@pytest.mark.integration
class TestCodonAdaptationIndex(object):
    def test_cai_basic(self, capsys):
        # Use the same CDS file as both query and reference. Each gene's
        # CAI should be in [0, 1].
        testargs = [
            "biokit", "codon_adaptation_index", SAMPLE, "-r", SAMPLE,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        assert out
        for line in out.split("\n"):
            fields = line.split("\t")
            assert len(fields) == 2
            if fields[1] != "None":
                value = float(fields[1])
                assert 0.0 <= value <= 1.0

    def test_cai_alias(self, capsys):
        testargs = ["biokit", "cai", SAMPLE, "-r", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        assert out

    def test_cai_verbose(self, capsys):
        testargs = ["biokit", "cai", SAMPLE, "-r", SAMPLE, "-v"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        for line in out.split("\n"):
            assert len(line.split("\t")) == 3

    def test_cai_json(self, capsys):
        testargs = [
            "biokit", "cai", SAMPLE, "-r", SAMPLE, "-f", "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        assert len(data) > 0
        for row in data:
            assert set(row.keys()) == {"gene_id", "CAI"}
