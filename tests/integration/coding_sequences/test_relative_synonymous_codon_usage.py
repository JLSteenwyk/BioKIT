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

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt1(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUAG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUAA\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "1",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt2(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nUUA\t2.1892\nAGA\t2.1667\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9355\nCGU\t1.8462\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAUU\t1.434\nAUG\t1.3913\nCAA\t1.3636\nAGG\t1.3333\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCGG\t1.2308\nCAU\t1.1667\nAAU\t1.0769\nACC\t1.0435\nGUG\t1.0169\nUAU\t1.0\nUAC\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nGUA\t0.7458\nUUC\t0.7442\nCUA\t0.6486\nCAG\t0.6364\nCGA\t0.6154\nAUA\t0.6087\nCUU\t0.5676\nCUG\t0.5676\nAUC\t0.566\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nUAG\t0.3333\nUGC\t0.3333\nGGC\t0.3158\nCGC\t0.3077\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUAA\t0.1667\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "2",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt3(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nACU\t2.1622\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9355\nGAA\t1.7714\nUGU\t1.6667\nAAA\t1.5\nAUU\t1.434\nAGG\t1.4118\nAUG\t1.3913\nCAA\t1.3636\nUAG\t1.3333\nACC\t1.2973\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nACA\t1.1892\nUUA\t1.1739\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nGUG\t1.0169\nUAU\t1.0\nUAC\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nGUC\t0.9492\nAAC\t0.9231\nCUA\t0.8649\nCCU\t0.8571\nCAC\t0.8333\nUUG\t0.8261\nGCA\t0.8148\nCUU\t0.7568\nCUG\t0.7568\nGAC\t0.7568\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nUAA\t0.6667\nCUC\t0.6486\nCAG\t0.6364\nAUA\t0.6087\nAUC\t0.566\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nUGC\t0.3333\nACG\t0.3243\nGGC\t0.3158\nCCC\t0.2857\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "3",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt4(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nUAG\t1.3333\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nUAA\t0.6667\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "4",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt5(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9277\nCGU\t1.8462\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAUU\t1.434\nAUG\t1.3913\nCAA\t1.3636\nUAG\t1.3333\nGUU\t1.2881\nUUU\t1.2558\nAGU\t1.253\nAGA\t1.253\nGAU\t1.2432\nCGG\t1.2308\nCAU\t1.1667\nAAU\t1.0769\nACC\t1.0435\nGUG\t1.0169\nUAU\t1.0\nUAC\t1.0\nUCC\t0.9639\nUCA\t0.9639\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nAGG\t0.7711\nGAC\t0.7568\nGUA\t0.7458\nUUC\t0.7442\nUAA\t0.6667\nCUA\t0.6486\nCAG\t0.6364\nCGA\t0.6154\nAUA\t0.6087\nCUU\t0.5676\nCUG\t0.5676\nAUC\t0.566\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4819\nCCG\t0.4286\nAGC\t0.3855\nUGC\t0.3333\nGGC\t0.3158\nCGC\t0.3077\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "5",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt6(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nCAA\t2.4\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nCAG\t1.12\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCUA\t0.6486\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nUAG\t0.32\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUAA\t0.16\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "6",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt9(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9277\nAUU\t1.9\nCGU\t1.8462\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.44\nCAA\t1.3636\nUAG\t1.3333\nGUU\t1.2881\nUUU\t1.2558\nAGU\t1.253\nAGA\t1.253\nGAU\t1.2432\nCGG\t1.2308\nCAU\t1.1667\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nAAG\t1.0\nUCC\t0.9639\nUCA\t0.9639\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nCCU\t0.8571\nAAU\t0.84\nCAC\t0.8333\nGCA\t0.8148\nAGG\t0.7711\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nAAC\t0.72\nUAA\t0.6667\nCUA\t0.6486\nCAG\t0.6364\nCGA\t0.6154\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nCUC\t0.4865\nUCG\t0.4819\nCCG\t0.4286\nAGC\t0.3855\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCGC\t0.3077\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "9",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt10(self, mocked_print):
        expected_result = "GGU\t2.6316\nUGU\t2.5\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nUAG\t1.3333\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nUAA\t0.6667\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nUGC\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "10",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt11(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUAG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUAA\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "11",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt12(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUCU\t2.029\nUUA\t2.0149\nGCU\t2.0\nUAG\t2.0\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nAAA\t1.5\nUUG\t1.4179\nAGG\t1.4118\nCAA\t1.3636\nAGU\t1.3188\nGUU\t1.2881\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nUCC\t1.0145\nUCA\t1.0145\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUAA\t1.0\nUGG\t1.0\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCUG\t0.7101\nCGG\t0.7059\nCAG\t0.6364\nCUA\t0.597\nGGA\t0.5263\nGGG\t0.5263\nCUU\t0.5224\nUCG\t0.5072\nAAG\t0.5\nCUC\t0.4478\nCCG\t0.4286\nAGC\t0.4058\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "12",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt13(self, mocked_print):
        expected_result = "GGU\t2.5424\nCCA\t2.4286\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9355\nCGU\t1.8462\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAUU\t1.434\nAUG\t1.3913\nCAA\t1.3636\nUAG\t1.3333\nAGA\t1.322\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCGG\t1.2308\nCAU\t1.1667\nAAU\t1.0769\nACC\t1.0435\nGUG\t1.0169\nUAU\t1.0\nUAC\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nAGG\t0.8136\nGAC\t0.7568\nGUA\t0.7458\nUUC\t0.7442\nUAA\t0.6667\nCUA\t0.6486\nCAG\t0.6364\nCGA\t0.6154\nAUA\t0.6087\nCUU\t0.5676\nCUG\t0.5676\nAUC\t0.566\nGGA\t0.5085\nGGG\t0.5085\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nUGC\t0.3333\nCGC\t0.3077\nGGC\t0.3051\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "13",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt14(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9277\nAUU\t1.9\nCGU\t1.8462\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nUAU\t1.4681\nUAC\t1.4681\nAAA\t1.44\nCAA\t1.3636\nGUU\t1.2881\nUUU\t1.2558\nAGU\t1.253\nAGA\t1.253\nGAU\t1.2432\nCGG\t1.2308\nCAU\t1.1667\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAG\t1.0\nAAG\t1.0\nUCC\t0.9639\nUCA\t0.9639\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nCCU\t0.8571\nAAU\t0.84\nCAC\t0.8333\nGCA\t0.8148\nAGG\t0.7711\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nAAC\t0.72\nCUA\t0.6486\nCAG\t0.6364\nCGA\t0.6154\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nCUC\t0.4865\nUCG\t0.4819\nCCG\t0.4286\nAGC\t0.3855\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCGC\t0.3077\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUAA\t0.0638\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "14",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt16(self, mocked_print):
        expected_result = "GGU\t2.6316\nUUA\t2.4868\nCCA\t2.4286\nAGA\t2.2941\nGCU\t2.0\nUAA\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nUUG\t1.75\nACU\t1.7391\nUGU\t1.6667\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCUA\t0.7368\nCGG\t0.7059\nCUU\t0.6447\nCUG\t0.6447\nCAG\t0.6364\nCUC\t0.5526\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUAG\t0.1842\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "16",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt21(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9277\nCGU\t1.8462\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.44\nAUU\t1.434\nAUG\t1.3913\nCAA\t1.3636\nUAG\t1.3333\nGUU\t1.2881\nUUU\t1.2558\nAGU\t1.253\nAGA\t1.253\nGAU\t1.2432\nCGG\t1.2308\nCAU\t1.1667\nACC\t1.0435\nGUG\t1.0169\nUAU\t1.0\nUAC\t1.0\nAAG\t1.0\nUCC\t0.9639\nUCA\t0.9639\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nCCU\t0.8571\nAAU\t0.84\nCAC\t0.8333\nGCA\t0.8148\nAGG\t0.7711\nGAC\t0.7568\nGUA\t0.7458\nUUC\t0.7442\nAAC\t0.72\nUAA\t0.6667\nCUA\t0.6486\nCAG\t0.6364\nCGA\t0.6154\nAUA\t0.6087\nCUU\t0.5676\nCUG\t0.5676\nAUC\t0.566\nGGA\t0.5263\nGGG\t0.5263\nCUC\t0.4865\nUCG\t0.4819\nCCG\t0.4286\nAGC\t0.3855\nUGC\t0.3333\nGGC\t0.3158\nCGC\t0.3077\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "21",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt22(self, mocked_print):
        expected_result = "UCA\t2.7273\nGGU\t2.6316\nUUA\t2.4868\nCCA\t2.4286\nAGA\t2.2941\nGCU\t2.0\nUCU\t1.9231\nAUU\t1.9\nGAA\t1.7714\nUUG\t1.75\nACU\t1.7391\nUGU\t1.6667\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nUUU\t1.2558\nAGU\t1.25\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUGG\t1.0\nGCC\t0.963\nUCC\t0.9615\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCUA\t0.7368\nCGG\t0.7059\nCUU\t0.6447\nCUG\t0.6447\nCAG\t0.6364\nCUC\t0.5526\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nUCG\t0.4808\nCCG\t0.4286\nAGC\t0.3846\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nUAA\t0.2727\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUAG\t0.1842\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "22",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt23(self, mocked_print):
        expected_result = "UUA\t3.6\nGGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUUG\t2.0213\nGCU\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCUA\t0.8511\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nCUU\t0.7447\nCUG\t0.7447\nUUC\t0.7442\nCGG\t0.7059\nCUC\t0.6383\nCAG\t0.6364\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nUAG\t0.2667\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUAA\t0.1333\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "23",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt24(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nAAA\t1.9286\nAUU\t1.9\nUCU\t1.8667\nCGU\t1.8462\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nCAA\t1.3636\nUAG\t1.3333\nGUU\t1.2881\nUUU\t1.2558\nGAU\t1.2432\nCGG\t1.2308\nAGU\t1.2133\nAGA\t1.2133\nCAU\t1.1667\nAAU\t1.0769\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nUCC\t0.9333\nUCA\t0.9333\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nUAA\t0.6667\nCUA\t0.6486\nAAG\t0.6429\nCAG\t0.6364\nCGA\t0.6154\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nCUC\t0.4865\nUCG\t0.4667\nCCG\t0.4286\nAGG\t0.4286\nAGC\t0.3733\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCGC\t0.3077\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "24",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt25(self, mocked_print):
        expected_result = "GGU\t3.2895\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nUAG\t1.3333\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nUAA\t0.6667\nGGA\t0.6579\nGGG\t0.6579\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nGGC\t0.3947\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "25",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt26(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nGCU\t2.2131\nUUA\t2.0149\nUAG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nAAA\t1.5\nUUG\t1.4179\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nGCC\t1.0656\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUAA\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nGCA\t0.9016\nCCU\t0.8571\nCAC\t0.8333\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCAG\t0.6364\nCUA\t0.597\nCUG\t0.5738\nGGA\t0.5263\nGGG\t0.5263\nCUU\t0.5224\nAAG\t0.5\nUCG\t0.4839\nCUC\t0.4478\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGCG\t0.2459\nGAG\t0.2286\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "26",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt27(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nCAA\t2.4\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nCAG\t1.12\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCUA\t0.6486\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nUAG\t0.32\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUAA\t0.16\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "27",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt28(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUAG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUAA\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "28",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt29(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUCU\t1.9355\nAUU\t1.9\nUAU\t1.8776\nUAC\t1.8776\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUAG\t0.1633\nUAA\t0.0816\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "29",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt30(self, mocked_print):
        expected_result = "GAA\t3.2632\nGGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUCU\t1.9355\nAUU\t1.9\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nGAG\t0.4211\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGCG\t0.2222\nUAG\t0.2105\nCGC\t0.1765\nUAA\t0.1053\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "30",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt31(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nAAA\t1.5\nAGG\t1.4118\nCAA\t1.3636\nUAG\t1.3333\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUCC\t0.9677\nUCA\t0.9677\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nUAA\t0.6667\nCUA\t0.6486\nCAG\t0.6364\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nAAG\t0.5\nCUC\t0.4865\nUCG\t0.4839\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "31",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt33(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nUUA\t2.1892\nGCU\t2.0\nUGG\t2.0\nAAA\t1.9286\nAUU\t1.9\nUCU\t1.8667\nCGU\t1.8462\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nUUG\t1.5405\nUAU\t1.4681\nUAC\t1.4681\nCAA\t1.3636\nGUU\t1.2881\nUUU\t1.2558\nGAU\t1.2432\nCGG\t1.2308\nAGU\t1.2133\nAGA\t1.2133\nCAU\t1.1667\nAAU\t1.0769\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAG\t1.0\nGCC\t0.963\nACA\t0.9565\nGUC\t0.9492\nUCC\t0.9333\nUCA\t0.9333\nAAC\t0.9231\nCCU\t0.8571\nCAC\t0.8333\nGCA\t0.8148\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCUA\t0.6486\nAAG\t0.6429\nCAG\t0.6364\nCGA\t0.6154\nCUU\t0.5676\nCUG\t0.5676\nGGA\t0.5263\nGGG\t0.5263\nCUC\t0.4865\nUCG\t0.4667\nCCG\t0.4286\nAGG\t0.4286\nAGC\t0.3733\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCGC\t0.3077\nCCC\t0.2857\nACG\t0.2609\nGAG\t0.2286\nGCG\t0.2222\nUAA\t0.0638\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "33",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt50(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nGCU\t2.2131\nUUA\t2.0149\nUAG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nAAA\t1.5\nUUG\t1.4179\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nGCC\t1.0656\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUAA\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nGCA\t0.9016\nCCU\t0.8571\nCAC\t0.8333\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCAG\t0.6364\nCUA\t0.597\nCUG\t0.5738\nGGA\t0.5263\nGGG\t0.5263\nCUU\t0.5224\nAAG\t0.5\nUCG\t0.4839\nCUC\t0.4478\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGCG\t0.2459\nGAG\t0.2286\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            "50",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt_custom(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nGCU\t2.2131\nUUA\t2.0149\nUAG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nAAA\t1.5\nUUG\t1.4179\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nGCC\t1.0656\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUAA\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nGCA\t0.9016\nCCU\t0.8571\nCAC\t0.8333\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCAG\t0.6364\nCUA\t0.597\nCUG\t0.5738\nGGA\t0.5263\nGGG\t0.5263\nCUU\t0.5224\nAAG\t0.5\nUCG\t0.4839\nCUC\t0.4478\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGCG\t0.2459\nGAG\t0.2286\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            f"{here.parent.parent.parent}/sample_files/CUG_ala_code.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_relative_synonymous_codon_usage_tt_custom_alias(self, mocked_print):
        expected_result = "GGU\t2.6316\nCCA\t2.4286\nAGA\t2.2941\nGCU\t2.2131\nUUA\t2.0149\nUAG\t2.0\nUCU\t1.9355\nAUU\t1.9\nGAA\t1.7714\nACU\t1.7391\nUGU\t1.6667\nAAA\t1.5\nUUG\t1.4179\nAGG\t1.4118\nCAA\t1.3636\nGUU\t1.2881\nAGU\t1.2581\nUUU\t1.2558\nGAU\t1.2432\nCAU\t1.1667\nAAU\t1.0769\nGCC\t1.0656\nCGU\t1.0588\nACC\t1.0435\nGUG\t1.0169\nAUG\t1.0\nUAU\t1.0\nUAC\t1.0\nUAA\t1.0\nUGG\t1.0\nUCC\t0.9677\nUCA\t0.9677\nACA\t0.9565\nGUC\t0.9492\nAAC\t0.9231\nGCA\t0.9016\nCCU\t0.8571\nCAC\t0.8333\nGAC\t0.7568\nAUC\t0.75\nGUA\t0.7458\nUUC\t0.7442\nCGG\t0.7059\nCAG\t0.6364\nCUA\t0.597\nCUG\t0.5738\nGGA\t0.5263\nGGG\t0.5263\nCUU\t0.5224\nAAG\t0.5\nUCG\t0.4839\nCUC\t0.4478\nCCG\t0.4286\nAGC\t0.3871\nCGA\t0.3529\nAUA\t0.35\nUGC\t0.3333\nGGC\t0.3158\nCCC\t0.2857\nACG\t0.2609\nGCG\t0.2459\nGAG\t0.2286\nCGC\t0.1765\nUGA\t0"

        testargs = [
            "biokit",
            "rscu",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.small.fna",
            "-tt",
            f"{here.parent.parent.parent}/sample_files/CUG_ala_code.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
