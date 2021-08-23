import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestRelativeSynonymousCodonUsage(object):
    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @pytest.mark.slow
    @patch("builtins.print")
    def test_relative_synonymous_codon_usage(self, mocked_print):
        expected_result = "AGA\t2.8477\nGGU\t1.8214\nUUA\t1.6781\nUUG\t1.6693\nCCA\t1.6278\nUCU\t1.5612\nGUU\t1.5466\nGCU\t1.4788\nGAA\t1.4019\nUAA\t1.3971\nAUU\t1.3842\nCAA\t1.3716\nACU\t1.3709\nGAU\t1.3034\nCAU\t1.2854\nAGG\t1.2769\nUCA\t1.2765\nUGU\t1.2441\nCCU\t1.2432\nACA\t1.2332\nAAU\t1.1952\nUUU\t1.1904\nGCA\t1.1875\nAAA\t1.1703\nUAU\t1.1366\nAUG\t1.0\nUGG\t1.0\nAGU\t0.9766\nUCC\t0.9386\nUGA\t0.929\nGGA\t0.9028\nGCC\t0.8833\nGUA\t0.8741\nUAC\t0.8634\nCUA\t0.853\nCGU\t0.848\nACC\t0.8439\nAUA\t0.84\nAAG\t0.8297\nUUC\t0.8096\nGUC\t0.8073\nAAC\t0.8048\nGGC\t0.7874\nCUU\t0.7781\nAUC\t0.7758\nGUG\t0.772\nUGC\t0.7559\nCAC\t0.7146\nGAC\t0.6966\nUAG\t0.6739\nCUG\t0.6719\nAGC\t0.6655\nCCC\t0.6323\nCAG\t0.6284\nGAG\t0.5981\nUCG\t0.5817\nACG\t0.552\nCCG\t0.4967\nGGG\t0.4885\nGCG\t0.4504\nCGA\t0.4208\nCGC\t0.358\nCUC\t0.3497\nCGG\t0.2485"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.fna",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
