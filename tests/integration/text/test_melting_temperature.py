import pytest

from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestMeltingTemperature(object):
    @patch("builtins.print")
    def test_melting_temperature(self, mocked_print):
        testargs = [
            "biokit",
            "melting_temperature",
            f"{here.parent.parent.parent}/sample_files/primers.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        # Should print one line per sequence (3 sequences in primers.fa)
        assert mocked_print.call_count == 3

    @patch("builtins.print")
    def test_melting_temperature_alias(self, mocked_print):
        testargs = [
            "biokit",
            "tm",
            f"{here.parent.parent.parent}/sample_files/primers.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        assert mocked_print.call_count == 3

    @patch("builtins.print")
    def test_melting_temperature_custom_na(self, mocked_print):
        testargs = [
            "biokit",
            "melting_temperature",
            f"{here.parent.parent.parent}/sample_files/primers.fa",
            "--na",
            "100",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        assert mocked_print.call_count == 3

    @patch("builtins.print")
    def test_melting_temperature_json(self, mocked_print):
        testargs = [
            "biokit",
            "melting_temperature",
            f"{here.parent.parent.parent}/sample_files/primers.fa",
            "-f",
            "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        import json
        output = mocked_print.call_args_list[0][0][0]
        data = json.loads(output)
        assert len(data) == 3
        assert all("Tm" in row for row in data)
