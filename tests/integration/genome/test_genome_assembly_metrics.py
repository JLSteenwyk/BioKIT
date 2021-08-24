import pytest
import textwrap

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestGenomeAssemblyMetrics(object):
    @patch("builtins.print")
    def test_genome_assembly_metrics_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @pytest.mark.slow
    @patch("builtins.print")
    def test_genome_assembly_metrics_simple(self, mocked_print):
        expected_result = textwrap.dedent(
            f"\
            12157105\tAssembly size\n\
            6\tL50\n\
            13\tL90\n\
            924431\tN50\n\
            439888\tN90\n\
            0.3815\tGC content\n\
            17\tNumber of scaffolds\n\
            17\tNumber of large scaffolds\n\
            12157105\tSum length of large scaffolds\n\
            0.3098\tFrequency of A\n\
            0.3087\tFrequency of T\n\
            0.1909\tFrequency of C\n\
            0.1906\tFrequency of G" # noqa
        )

        testargs = [
            "biokit",
            "genome_assembly_metrics",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @pytest.mark.slow
    @patch("builtins.print")
    def test_genome_assembly_metrics_simple_alias(self, mocked_print):
        expected_result = textwrap.dedent(
            f"\
            12157105\tAssembly size\n\
            6\tL50\n\
            13\tL90\n\
            924431\tN50\n\
            439888\tN90\n\
            0.3815\tGC content\n\
            17\tNumber of scaffolds\n\
            17\tNumber of large scaffolds\n\
            12157105\tSum length of large scaffolds\n\
            1531933\tLongest scaffold\n\
            0.3098\tFrequency of A\n\
            0.3087\tFrequency of T\n\
            0.1909\tFrequency of C\n\
            0.1906\tFrequency of G" # noqa
        )

        testargs = [
            "biokit",
            "assembly_metrics",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
