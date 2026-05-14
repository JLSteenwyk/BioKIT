import json
import os

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
class TestEffectiveNumberOfCodons(object):
    def test_enc_basic(self, capsys):
        testargs = ["biokit", "effective_number_of_codons", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        lines = out.split("\n")
        # Each line has 2 columns: gene_id and ENC
        for line in lines:
            fields = line.split("\t")
            assert len(fields) == 2
            # ENC must be a number in [20, 61]
            value = float(fields[1])
            assert 20.0 <= value <= 61.0

    def test_enc_alias(self, capsys):
        testargs = ["biokit", "enc", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        assert out  # non-empty

    def test_enc_verbose(self, capsys):
        testargs = ["biokit", "enc", SAMPLE, "-v"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        # Each line should now have 3 columns
        for line in out.split("\n"):
            assert len(line.split("\t")) == 3

    def test_enc_json(self, capsys):
        testargs = ["biokit", "enc", SAMPLE, "-f", "json"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        assert len(data) > 0
        for row in data:
            assert set(row.keys()) == {"gene_id", "ENC"}

    def test_enc_plot(self, tmp_path):
        out_path = str(tmp_path / "enc.png")
        testargs = [
            "biokit", "enc", SAMPLE, "--plot", "-o", out_path,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert os.path.exists(out_path)
        assert os.path.getsize(out_path) > 0
