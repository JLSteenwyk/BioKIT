import json

import pytest
from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestHomopolymerRuns(object):
    def test_homopolymer_runs_basic(self, capsys):
        testargs = [
            "biokit",
            "homopolymer_runs",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert lines[0] == "id\tlength\tbase\tposition"
        # at least one sequence row plus header
        assert len(lines) > 1

    def test_homopolymer_runs_alias(self, capsys):
        testargs = [
            "biokit",
            "homopolymer",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        assert out.startswith("id\tlength\tbase\tposition")

    def test_homopolymer_runs_json(self, capsys):
        testargs = [
            "biokit",
            "homopolymer_runs",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "-f",
            "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        data = json.loads(out)
        assert len(data) > 0
        for row in data:
            assert set(row.keys()) == {"id", "length", "base", "position"}

    def test_homopolymer_runs_per_base(self, capsys):
        testargs = [
            "biokit",
            "homopolymer_runs",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "--per-base",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert lines[0] == "id\tA\tC\tG\tT"
