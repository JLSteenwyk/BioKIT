import json
import os

import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna"


@pytest.mark.integration
class TestNeutralityPlot(object):
    def test_neutrality_plot_summary(self, capsys):
        testargs = ["biokit", "neutrality_plot", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        # 4 summary lines
        lines = out.strip().split("\n")
        labels = [line.split("\t")[1] for line in lines]
        assert "slope" in labels
        assert "intercept" in labels
        assert "r_squared" in labels
        assert "n" in labels

    def test_neutrality_plot_verbose(self, capsys):
        testargs = ["biokit", "neutrality_plot", SAMPLE, "-v"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert "id\tGC12\tGC3" in out
        assert "slope" in out

    def test_neutrality_plot_json(self, capsys):
        testargs = ["biokit", "neutrality_plot", SAMPLE, "-f", "json"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        assert "regression" in data
        assert "n" in data["regression"]

    def test_neutrality_plot_plot(self, tmp_path):
        out_path = str(tmp_path / "plot.png")
        testargs = [
            "biokit", "neutrality_plot", SAMPLE,
            "--plot", "-o", out_path,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert os.path.exists(out_path)
        assert os.path.getsize(out_path) > 0
