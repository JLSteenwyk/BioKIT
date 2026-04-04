import pytest

from mock import patch
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestGCContentFourFoldDegenerateSites(object):
    @patch("builtins.print")
    def test_gc4_invalid_input(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_gc4_simple(self, mocked_print):
        testargs = [
            "biokit",
            "gc_content_four_fold_degenerate_sites",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        result = mocked_print.mock_calls[0][1][0]
        # Should return a numeric string
        float(result)

    @patch("builtins.print")
    def test_gc4_verbose(self, mocked_print):
        testargs = [
            "biokit",
            "gc_content_four_fold_degenerate_sites",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-v",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        # Should have 3 calls (one per sequence in the small file)
        assert len(mocked_print.mock_calls) == 3
        # Each call should have seq_id\tvalue format
        for c in mocked_print.mock_calls:
            line = c[1][0]
            parts = line.split("\t")
            assert len(parts) == 2
            float(parts[1])

    @patch("builtins.print")
    def test_gc4_alias(self, mocked_print):
        testargs = [
            "biokit",
            "gc4",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-v",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()

        assert len(mocked_print.mock_calls) == 3
