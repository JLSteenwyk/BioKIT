import pytest

from mock import patch, call  # noqa
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestProteinCharge(object):
    @patch("builtins.print")
    def test_protein_charge(self, mocked_print):
        testargs = [
            "biokit",
            "protein_charge",
            f"{here.parent.parent.parent}/sample_files/protein.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        output = mocked_print.call_args_list[0][0][0]
        assert "insulin_A" in output

    @patch("builtins.print")
    def test_protein_charge_alias(self, mocked_print):
        testargs = [
            "biokit",
            "prot_charge",
            f"{here.parent.parent.parent}/sample_files/protein.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        output = mocked_print.call_args_list[0][0][0]
        assert "insulin_A" in output

    @patch("builtins.print")
    def test_protein_charge_custom_pH(self, mocked_print):
        testargs = [
            "biokit",
            "protein_charge",
            f"{here.parent.parent.parent}/sample_files/protein.fa",
            "-p",
            "4.0",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        output = mocked_print.call_args_list[0][0][0]
        assert "insulin_A" in output

    @patch("builtins.print")
    def test_protein_charge_json_format(self, mocked_print):
        testargs = [
            "biokit",
            "protein_charge",
            f"{here.parent.parent.parent}/sample_files/protein.fa",
            "-f",
            "json",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        import json
        output = mocked_print.call_args_list[0][0][0]
        data = json.loads(output)
        assert len(data) == 3
        assert any(row["id"] == "insulin_A" for row in data)
