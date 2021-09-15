import pytest

from mock import patch, call
from pathlib import Path
import sys

from biokit.biokit import Biokit

here = Path(__file__)


@pytest.mark.integration
class TestGeneWiseRelativeSynonymousCodonUsage(object):
    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Biokit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1543\t1.2233\t0.4166\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.166\t1.2411\t0.4201\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1321\t1.2233\t0.3879\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1782\t1.2411\t0.4411\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2074\t1.2749\t0.4157\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1584\t1.2411\t0.4245\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1528\t1.2411\t0.419\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1722\t1.2411\t0.4286\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1532\t1.2411\t0.4178\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.167\t1.2411\t0.4159\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1703\t1.2411\t0.4139\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1522\t1.2411\t0.4257"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt1(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1543\t1.2233\t0.4166\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.166\t1.2411\t0.4201\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1321\t1.2233\t0.3879\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1782\t1.2411\t0.4411\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2074\t1.2749\t0.4157\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1584\t1.2411\t0.4245\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1528\t1.2411\t0.419\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1722\t1.2411\t0.4286\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1532\t1.2411\t0.4178\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.167\t1.2411\t0.4159\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1703\t1.2411\t0.4139\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1522\t1.2411\t0.4257"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "1",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt2(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1751\t1.2411\t0.424\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1774\t1.2411\t0.4238\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1396\t1.2411\t0.3906\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1954\t1.2587\t0.4469\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2245\t1.3265\t0.4233\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1806\t1.2411\t0.4346\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1724\t1.2411\t0.4257\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.188\t1.2411\t0.4351\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1757\t1.2411\t0.4283\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.188\t1.2587\t0.426\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.187\t1.2411\t0.4226\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1665\t1.2411\t0.4291"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "2",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt3(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1479\t1.2233\t0.394\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1449\t1.2233\t0.3953\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1216\t1.2233\t0.362\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1606\t1.2411\t0.4179\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1726\t1.2411\t0.3888\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1449\t1.2233\t0.3983\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1502\t1.2411\t0.4032\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1508\t1.2233\t0.4085\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1506\t1.2233\t0.4003\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1546\t1.2411\t0.3951\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1591\t1.2411\t0.3938\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.134\t1.2233\t0.4009"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "3",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt4(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1688\t1.2411\t0.428\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1732\t1.2411\t0.4257\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1349\t1.2411\t0.3906\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1846\t1.2411\t0.4462\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2126\t1.2952\t0.42\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1692\t1.2411\t0.4327\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1655\t1.2411\t0.4291\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1817\t1.2411\t0.4361\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1674\t1.2411\t0.429\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1776\t1.2411\t0.424\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1815\t1.2411\t0.4224\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1585\t1.2411\t0.4307"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "4",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt5(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1468\t1.2411\t0.3763\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1477\t1.2411\t0.3749\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1253\t1.2411\t0.3636\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.157\t1.2411\t0.3817\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1989\t1.2749\t0.3887\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1496\t1.2411\t0.3868\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.141\t1.2411\t0.3778\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1562\t1.2411\t0.3854\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1489\t1.2411\t0.3831\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1626\t1.2411\t0.3848\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.156\t1.2411\t0.3711\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1365\t1.2411\t0.3804"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "5",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt6(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.2\t1.2233\t0.4932\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1959\t1.2411\t0.4703\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1704\t1.2233\t0.4683\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.212\t1.2411\t0.4952\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2459\t1.2749\t0.4784\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.205\t1.2411\t0.4982\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.2022\t1.2411\t0.5036\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.2189\t1.2411\t0.4998\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1855\t1.2411\t0.4706\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.2091\t1.2411\t0.487\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.2087\t1.2411\t0.4787\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1884\t1.2411\t0.4876"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "6",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt9(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1452\t1.2463\t0.3785\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1447\t1.2463\t0.3772\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1309\t1.2463\t0.3609\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1575\t1.2587\t0.3832\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2017\t1.2946\t0.3885\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1476\t1.2463\t0.389\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.142\t1.2463\t0.3785\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.157\t1.2587\t0.3863\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1498\t1.2587\t0.3843\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1637\t1.2587\t0.3851\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1566\t1.2587\t0.3727\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.135\t1.2463\t0.3825"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "9",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt10(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1629\t1.2411\t0.4194\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.173\t1.2411\t0.4243\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1371\t1.2411\t0.3892\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1838\t1.2411\t0.4442\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2117\t1.2952\t0.4168\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1634\t1.2411\t0.4266\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1591\t1.2411\t0.4223\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1756\t1.2411\t0.4311\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1584\t1.2411\t0.4203\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1708\t1.2411\t0.4169\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1759\t1.2411\t0.4169\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1586\t1.2411\t0.4289"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "10",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt11(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1543\t1.2233\t0.4166\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.166\t1.2411\t0.4201\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1321\t1.2233\t0.3879\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1782\t1.2411\t0.4411\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2074\t1.2749\t0.4157\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1584\t1.2411\t0.4245\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1528\t1.2411\t0.419\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1722\t1.2411\t0.4286\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1532\t1.2411\t0.4178\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.167\t1.2411\t0.4159\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1703\t1.2411\t0.4139\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1522\t1.2411\t0.4257"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "11",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt12(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1484\t1.2233\t0.4063\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1612\t1.2411\t0.4111\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1289\t1.2233\t0.3797\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1728\t1.2411\t0.4325\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1992\t1.2749\t0.4044\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1516\t1.2411\t0.4137\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1475\t1.2411\t0.4115\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1677\t1.2411\t0.419\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1506\t1.2411\t0.4089\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1611\t1.2411\t0.4071\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1646\t1.2411\t0.4053\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1477\t1.2411\t0.4171"
        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "12",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt13(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.149\t1.2411\t0.38\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1508\t1.2411\t0.3782\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1251\t1.2411\t0.3658\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1621\t1.2411\t0.3873\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1971\t1.3265\t0.3862\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1548\t1.2411\t0.3903\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1439\t1.2411\t0.3795\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.158\t1.2411\t0.3879\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1492\t1.2411\t0.3858\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1623\t1.2411\t0.3861\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1588\t1.2411\t0.3753\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1398\t1.2411\t0.3835"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "13",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt14(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1643\t1.2587\t0.3864\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1623\t1.2587\t0.3865\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1515\t1.2749\t0.3696\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1756\t1.2749\t0.392\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2114\t1.3265\t0.3912\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1649\t1.2587\t0.3951\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1578\t1.2587\t0.388\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1715\t1.2749\t0.3911\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1684\t1.2587\t0.3935\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1778\t1.2749\t0.3913\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1703\t1.2587\t0.3803\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.153\t1.2587\t0.3905"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "14",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt16(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.181\t1.2233\t0.4423\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1918\t1.2411\t0.4451\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1546\t1.2233\t0.4087\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.2023\t1.2411\t0.4626\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2363\t1.2749\t0.4442\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1853\t1.2411\t0.4504\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1765\t1.2411\t0.4382\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1951\t1.2411\t0.4528\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1736\t1.2411\t0.438\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.191\t1.2411\t0.4389\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1935\t1.2411\t0.4366\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1752\t1.2411\t0.4497"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "16",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt21(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1455\t1.2463\t0.3766\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1443\t1.2463\t0.3753\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1302\t1.2463\t0.3596\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1571\t1.2587\t0.3815\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2\t1.2946\t0.3868\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1477\t1.2463\t0.3872\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.142\t1.2463\t0.3765\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1559\t1.2587\t0.3848\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1497\t1.2587\t0.3823\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1628\t1.2587\t0.3838\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1561\t1.2587\t0.3708\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1342\t1.2463\t0.3805"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "21",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt22(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.222\t1.2233\t0.5161\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.2511\t1.2411\t0.5485\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1965\t1.2233\t0.4882\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.2394\t1.2411\t0.525\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.265\t1.2749\t0.491\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.2187\t1.2411\t0.5099\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.2093\t1.2411\t0.4965\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.2326\t1.2411\t0.5159\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.2193\t1.2411\t0.5214\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.2258\t1.2411\t0.5005\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.2274\t1.2411\t0.4974\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.2249\t1.2411\t0.5354"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "22",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt23(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.2586\t1.2233\t0.6873\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.2623\t1.2411\t0.6688\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.2078\t1.2233\t0.5957\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.2717\t1.2411\t0.6768\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.321\t1.2749\t0.6996\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.2636\t1.2411\t0.692\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.2313\t1.2411\t0.6184\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.2649\t1.2411\t0.6736\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.2357\t1.2411\t0.6451\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.259\t1.2411\t0.6581\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.2609\t1.2411\t0.656\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.2435\t1.2411\t0.6674"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "23",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt24(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.165\t1.2233\t0.3934\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1596\t1.2233\t0.3882\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1533\t1.2233\t0.3768\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1721\t1.2587\t0.3937\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2147\t1.3397\t0.3981\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1645\t1.2233\t0.4006\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1588\t1.2233\t0.3941\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1738\t1.2587\t0.3984\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1654\t1.2587\t0.3972\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1824\t1.2587\t0.3974\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.181\t1.2587\t0.3891\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1475\t1.2233\t0.3938"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "24",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt25(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1643\t1.2411\t0.4235\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1743\t1.2411\t0.4255\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1377\t1.2411\t0.3903\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1889\t1.2411\t0.4505\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2228\t1.2749\t0.4338\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1653\t1.2411\t0.43\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1638\t1.2411\t0.4296\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1844\t1.2411\t0.4395\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1655\t1.2411\t0.4276\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1799\t1.2411\t0.4278\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1813\t1.2411\t0.4212\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1606\t1.2411\t0.4323"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "25",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt26(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1473\t1.2233\t0.4041\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1589\t1.2411\t0.4085\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.128\t1.2233\t0.377\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1719\t1.2411\t0.4308\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1997\t1.2749\t0.4033\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1517\t1.2411\t0.4119\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1468\t1.2411\t0.4096\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1662\t1.2411\t0.4172\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1493\t1.2411\t0.4069\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1605\t1.2411\t0.4057\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1632\t1.2411\t0.4037\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1455\t1.2411\t0.4151"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "26",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt27(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.2\t1.2233\t0.4932\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1959\t1.2411\t0.4703\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1704\t1.2233\t0.4683\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.212\t1.2411\t0.4952\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2459\t1.2749\t0.4784\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.205\t1.2411\t0.4982\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.2022\t1.2411\t0.5036\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.2189\t1.2411\t0.4998\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1855\t1.2411\t0.4706\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.2091\t1.2411\t0.487\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.2087\t1.2411\t0.4787\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1884\t1.2411\t0.4876"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "27",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt28(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1543\t1.2233\t0.4166\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.166\t1.2411\t0.4201\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1321\t1.2233\t0.3879\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1782\t1.2411\t0.4411\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2074\t1.2749\t0.4157\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1584\t1.2411\t0.4245\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1528\t1.2411\t0.419\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1722\t1.2411\t0.4286\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1532\t1.2411\t0.4178\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.167\t1.2411\t0.4159\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1703\t1.2411\t0.4139\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1522\t1.2411\t0.4257"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "28",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt29(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1928\t1.2411\t0.4513\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.2015\t1.2587\t0.455\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1729\t1.2411\t0.4287\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.2144\t1.2587\t0.4747\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2268\t1.3265\t0.432\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1926\t1.2587\t0.4542\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1839\t1.2411\t0.4543\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.2013\t1.2587\t0.4523\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1906\t1.2411\t0.4554\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1952\t1.2587\t0.4435\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1979\t1.2411\t0.4427\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1876\t1.2587\t0.4591"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "29",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt30(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.2261\t1.2233\t0.5256\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.2285\t1.2411\t0.5159\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.2633\t1.2233\t0.5909\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.2625\t1.2411\t0.5581\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.3005\t1.2749\t0.5519\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.2359\t1.2411\t0.5469\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.2391\t1.2411\t0.5448\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.2569\t1.2411\t0.5531\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.2309\t1.2411\t0.5401\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.2727\t1.2411\t0.5754\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.2404\t1.2411\t0.5208\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.2196\t1.2411\t0.5323"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "30",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt31(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1688\t1.2411\t0.428\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1732\t1.2411\t0.4257\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.1349\t1.2411\t0.3906\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1846\t1.2411\t0.4462\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2126\t1.2952\t0.42\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1692\t1.2411\t0.4327\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1655\t1.2411\t0.4291\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1817\t1.2411\t0.4361\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1674\t1.2411\t0.429\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1776\t1.2411\t0.424\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1815\t1.2411\t0.4224\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1585\t1.2411\t0.4307"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "31",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt33(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1841\t1.2587\t0.4\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1773\t1.2587\t0.3966\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.174\t1.2749\t0.3839\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1902\t1.2749\t0.4016\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.2244\t1.3627\t0.4004\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1817\t1.2587\t0.4058\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1745\t1.2587\t0.4026\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1883\t1.2749\t0.4025\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.184\t1.2587\t0.4054\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1964\t1.2749\t0.4028\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1946\t1.2587\t0.3955\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1654\t1.2587\t0.401"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "33",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt50(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1473\t1.2233\t0.4041\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1589\t1.2411\t0.4085\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.128\t1.2233\t0.377\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1719\t1.2411\t0.4308\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1997\t1.2749\t0.4033\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1517\t1.2411\t0.4119\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1468\t1.2411\t0.4096\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1662\t1.2411\t0.4172\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1493\t1.2411\t0.4069\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1605\t1.2411\t0.4057\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1632\t1.2411\t0.4037\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1455\t1.2411\t0.4151"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            "50",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt_custom(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1473\t1.2233\t0.4041\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1589\t1.2411\t0.4085\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.128\t1.2233\t0.377\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1719\t1.2411\t0.4308\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1997\t1.2749\t0.4033\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1517\t1.2411\t0.4119\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1468\t1.2411\t0.4096\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1662\t1.2411\t0.4172\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1493\t1.2411\t0.4069\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1605\t1.2411\t0.4057\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1632\t1.2411\t0.4037\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1455\t1.2411\t0.4151"

        testargs = [
            "biokit",
            "gene_wise_relative_synonymous_codon_usage",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            f"{here.parent.parent.parent}/sample_files/CUG_ala_code.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt_custom_alias0(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1473\t1.2233\t0.4041\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1589\t1.2411\t0.4085\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.128\t1.2233\t0.377\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1719\t1.2411\t0.4308\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1997\t1.2749\t0.4033\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1517\t1.2411\t0.4119\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1468\t1.2411\t0.4096\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1662\t1.2411\t0.4172\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1493\t1.2411\t0.4069\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1605\t1.2411\t0.4057\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1632\t1.2411\t0.4037\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1455\t1.2411\t0.4151"

        testargs = [
            "biokit",
            "gene_wise_rscu",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            f"{here.parent.parent.parent}/sample_files/CUG_ala_code.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gene_wise_relative_synonymous_codon_usage_tt_custom_alias1(self, mocked_print):
        expected_result = "lcl|NC_001134.8_cds_NP_009465.2_119\t1.1473\t1.2233\t0.4041\nlcl|NC_001134.8_cds_NP_009698.3_345\t1.1589\t1.2411\t0.4085\nlcl|NC_001136.10_cds_NP_010434.1_1057\t1.128\t1.2233\t0.377\nlcl|NC_001136.10_cds_NP_010745.3_1362\t1.1719\t1.2411\t0.4308\nlcl|NC_001139.9_cds_NP_011320.3_1926\t1.1997\t1.2749\t0.4033\nlcl|NC_001140.6_cds_NP_011967.1_2560\t1.1517\t1.2411\t0.4119\nlcl|NC_001143.9_cds_NP_012980.1_3530\t1.1468\t1.2411\t0.4096\nlcl|NC_001144.5_cds_NP_013060.1_3608\t1.1662\t1.2411\t0.4172\nlcl|NC_001144.5_cds_NP_013188.1_3730\t1.1493\t1.2411\t0.4069\nlcl|NC_001144.5_cds_NP_013207.1_3749\t1.1605\t1.2411\t0.4057\nlcl|NC_001144.5_cds_NP_013559.1_4092\t1.1632\t1.2411\t0.4037\nlcl|NC_001147.6_cds_NP_014560.1_5058\t1.1455\t1.2411\t0.4151"

        testargs = [
            "biokit",
            "gw_rscu",
            f"{here.parent.parent.parent}/sample_files/GCF_000146045.2_R64_cds_from_genomic.long_sequences.fa.long_seqs.fa",
            "-tt",
            f"{here.parent.parent.parent}/sample_files/CUG_ala_code.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Biokit()
        assert mocked_print.mock_calls == [call(expected_result)]
