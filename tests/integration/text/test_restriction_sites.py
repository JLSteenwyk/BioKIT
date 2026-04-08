import pytest

from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestRestrictionSites(object):
    @patch("builtins.print")
    def test_restriction_sites(self, mocked_print):
        testargs = [
            "biokit",
            "restriction_sites",
            f"{here.parent.parent.parent}/sample_files/restriction_test.fa",
            "-e",
            "EcoRI",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        calls = [str(c) for c in mocked_print.call_args_list]
        output = "\n".join(calls)
        assert "EcoRI" in output

    @patch("builtins.print")
    def test_restriction_sites_alias(self, mocked_print):
        testargs = [
            "biokit",
            "re_sites",
            f"{here.parent.parent.parent}/sample_files/restriction_test.fa",
            "-e",
            "EcoRI,BamHI",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        calls = [str(c) for c in mocked_print.call_args_list]
        output = "\n".join(calls)
        assert "EcoRI" in output
        assert "BamHI" in output

    @patch("builtins.print")
    def test_restriction_sites_multiple_e_flags(self, mocked_print):
        testargs = [
            "biokit",
            "restriction_sites",
            f"{here.parent.parent.parent}/sample_files/restriction_test.fa",
            "-e",
            "EcoRI",
            "-e",
            "BamHI",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        calls = [str(c) for c in mocked_print.call_args_list]
        output = "\n".join(calls)
        assert "EcoRI" in output
        assert "BamHI" in output
