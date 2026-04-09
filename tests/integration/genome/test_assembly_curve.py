import json
import os

import pytest

from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestAssemblyCurve(object):
    @patch("builtins.print")
    def test_assembly_curve(self, mocked_print):
        testargs = [
            "biokit",
            "assembly_curve",
            f"{here.parent.parent.parent}/sample_files/test_assembly.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        # Header + 5 contigs = 6 print calls
        assert mocked_print.call_count == 6

    @patch("builtins.print")
    def test_assembly_curve_alias(self, mocked_print):
        testargs = [
            "biokit",
            "asm_curve",
            f"{here.parent.parent.parent}/sample_files/test_assembly.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        assert mocked_print.call_count == 6

    @patch("builtins.print")
    def test_assembly_curve_json(self, mocked_print):
        testargs = [
            "biokit",
            "assembly_curve",
            f"{here.parent.parent.parent}/sample_files/test_assembly.fa",
            "-f",
            "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        output = mocked_print.call_args_list[0][0][0]
        data = json.loads(output)
        assert len(data) == 5
        assert all("rank" in row for row in data)
        assert all("cumulative_length" in row for row in data)

    @patch("builtins.print")
    def test_assembly_curve_plot(self, mocked_print, tmp_path):
        plot_path = str(tmp_path / "test_curve.png")
        testargs = [
            "biokit",
            "assembly_curve",
            f"{here.parent.parent.parent}/sample_files/test_assembly.fa",
            "--plot",
            "-o",
            plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0
