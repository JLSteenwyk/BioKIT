import json

import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa"


@pytest.mark.integration
class TestTransitionTransversionRatio(object):
    def test_ti_tv_basic(self, capsys):
        testargs = ["biokit", "transition_transversion_ratio", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        fields = out.split("\t")
        assert len(fields) == 3
        # transitions and transversions parse as integers
        int(fields[0])
        int(fields[1])

    def test_ti_tv_alias(self, capsys):
        testargs = ["biokit", "ti_tv", SAMPLE]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        assert len(out.split("\t")) == 3

    def test_ti_tv_verbose(self, capsys):
        testargs = ["biokit", "ti_tv", SAMPLE, "-v"]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out.strip()
        labels = {line.split("\t")[1] for line in out.split("\n")}
        assert labels.issubset(
            {"transition", "transversion", "constant", "gap"}
        )

    def test_ti_tv_json(self, capsys):
        testargs = [
            "biokit", "ti_tv", SAMPLE, "-f", "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        assert set(data.keys()) == {"transitions", "transversions", "ratio"}

    def test_ti_tv_json_verbose(self, capsys):
        testargs = [
            "biokit", "ti_tv", SAMPLE, "-v", "-f", "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        data = json.loads(capsys.readouterr().out)
        for row in data:
            assert set(row.keys()) == {"position", "classification"}
