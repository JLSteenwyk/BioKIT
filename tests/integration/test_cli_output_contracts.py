import json
import sys
from pathlib import Path

import pytest
from mock import patch

from biokit.biokit import Biokit

here = Path(__file__)
sample_files = here.parent.parent / "sample_files"

GENOME_ASSEMBLY_LABEL_TO_KEY = {
    "Assembly size": "assembly_size",
    "L50": "l50",
    "L90": "l90",
    "N50": "n50",
    "N90": "n90",
    "GC content": "gc_content",
    "Number of scaffolds": "number_of_scaffolds",
    "Number of large scaffolds": "number_of_large_scaffolds",
    "Sum length of large scaffolds": "sum_length_of_large_scaffolds",
    "Longest scaffold": "longest_scaffold",
    "Frequency of A": "frequency_of_a",
    "Frequency of T": "frequency_of_t",
    "Frequency of C": "frequency_of_c",
    "Frequency of G": "frequency_of_g",
}


def _run_cli(args: list[str]) -> str:
    with patch("builtins.print") as mocked_print:
        with patch.object(sys, "argv", args):
            Biokit()

    return "\n".join(str(call.args[0]) for call in mocked_print.call_args_list)


def _parse_scalar(value: str):
    if value == "null":
        return None
    if value == "true":
        return True
    if value == "false":
        return False
    if value.startswith('"') and value.endswith('"'):
        return json.loads(value)
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _parse_yaml_simple(text: str):
    lines = [line.rstrip("\n") for line in text.splitlines() if line.strip()]
    if not lines:
        return {}

    if lines[0].startswith("- "):
        parsed = []
        current = None
        for line in lines:
            if line.startswith("- "):
                if current is not None:
                    parsed.append(current)
                current = {}
                head = line[2:]
                key, value = head.split(": ", 1)
                current[key] = _parse_scalar(value)
            else:
                key, value = line.strip().split(": ", 1)
                current[key] = _parse_scalar(value)
        if current is not None:
            parsed.append(current)
        return parsed

    parsed = {}
    for line in lines:
        key, value = line.split(": ", 1)
        parsed[key] = _parse_scalar(value)
    return parsed


def _run_and_parse(
    command: str,
    target_file: Path,
    output_format: str,
    extra_args: list[str] | None = None,
):
    extra_args = extra_args or []
    output = _run_cli(
        ["biokit", command, str(target_file), *extra_args, "--format", output_format]
    )

    if output_format == "json":
        return json.loads(output)
    if output_format == "yaml":
        return _parse_yaml_simple(output)
    return output


def _rows_to_map(rows, key_field: str):
    mapped = {}
    for row in rows:
        key = row[key_field]
        mapped[key] = {k: v for k, v in row.items() if k != key_field}
    return mapped


@pytest.mark.integration
class TestCliOutputContracts(object):
    def test_sequence_length_contract(self):
        target = sample_files / "EOG091N44MS.fa"

        tsv_output = _run_and_parse("sequence_length", target, "tsv")
        json_output = _run_and_parse("sequence_length", target, "json")
        yaml_output = _run_and_parse("sequence_length", target, "yaml")

        tsv_rows = []
        for row in tsv_output.splitlines():
            entry_id, length = row.split("\t")
            tsv_rows.append({"id": entry_id, "length": int(length)})
        tsv_rows.sort(key=lambda row: row["id"])

        assert tsv_rows == json_output
        assert tsv_rows == yaml_output
        assert list(json_output[0].keys()) == ["id", "length"]
        assert list(yaml_output[0].keys()) == ["id", "length"]

    def test_character_frequency_contract(self):
        target = sample_files / "simple.fa"

        tsv_output = _run_and_parse("character_frequency", target, "tsv")
        json_output = _run_and_parse("character_frequency", target, "json")
        yaml_output = _run_and_parse("character_frequency", target, "yaml")

        tsv_map = {}
        for row in tsv_output.splitlines():
            char, freq = row.split("\t")
            tsv_map[char] = float(freq)

        json_map = {row["character"]: row["frequency"] for row in json_output}
        yaml_map = {row["character"]: row["frequency"] for row in yaml_output}

        assert tsv_map == json_map
        assert tsv_map == yaml_map
        assert list(json_output[0].keys()) == ["character", "frequency"]
        assert list(yaml_output[0].keys()) == ["character", "frequency"]
        assert [row["character"] for row in json_output] == sorted(json_map.keys())

    def test_genome_assembly_metrics_contract(self):
        target = sample_files / "simple.fa"

        tsv_output = _run_and_parse("genome_assembly_metrics", target, "tsv")
        json_output = _run_and_parse("genome_assembly_metrics", target, "json")
        yaml_output = _run_and_parse("genome_assembly_metrics", target, "yaml")

        tsv_map = {}
        for row in tsv_output.splitlines():
            value, label = row.split("\t")
            key = GENOME_ASSEMBLY_LABEL_TO_KEY[label]
            tsv_map[key] = _parse_scalar(value)

        expected_keys = [
            "assembly_size",
            "l50",
            "l90",
            "n50",
            "n90",
            "gc_content",
            "number_of_scaffolds",
            "number_of_large_scaffolds",
            "sum_length_of_large_scaffolds",
            "longest_scaffold",
            "frequency_of_a",
            "frequency_of_t",
            "frequency_of_c",
            "frequency_of_g",
        ]

        assert tsv_map == json_output
        assert tsv_map == yaml_output
        assert list(json_output.keys()) == expected_keys
        assert list(yaml_output.keys()) == expected_keys

    def test_gc_content_verbose_contract(self):
        target = sample_files / "simple.fa"

        tsv_output = _run_and_parse("gc_content", target, "tsv", ["-v"])
        json_output = _run_and_parse("gc_content", target, "json", ["-v"])
        yaml_output = _run_and_parse("gc_content", target, "yaml", ["-v"])

        tsv_map = {}
        for row in tsv_output.splitlines():
            entry, gc_content = row.split("\t")
            tsv_map[entry] = float(gc_content)

        json_map = {str(row["entry"]): row["gc_content"] for row in json_output}
        yaml_map = {str(row["entry"]): row["gc_content"] for row in yaml_output}

        assert tsv_map == json_map
        assert tsv_map == yaml_map
        assert list(json_output[0].keys()) == ["entry", "gc_content"]
        assert list(yaml_output[0].keys()) == ["entry", "gc_content"]

    def test_fastq_read_lengths_verbose_contract(self):
        target = sample_files / "DRR284700_1_subset.fastq"

        tsv_output = _run_and_parse("fastq_read_lengths", target, "tsv", ["-v"])
        json_output = _run_and_parse("fastq_read_lengths", target, "json", ["-v"])
        yaml_output = _run_and_parse("fastq_read_lengths", target, "yaml", ["-v"])

        tsv_lengths = [int(line) for line in tsv_output.splitlines()]
        json_lengths = [row["length"] for row in json_output]
        yaml_lengths = [row["length"] for row in yaml_output]

        assert tsv_lengths == json_lengths
        assert tsv_lengths == yaml_lengths
        assert list(json_output[0].keys()) == ["read_index", "length"]
        assert list(yaml_output[0].keys()) == ["read_index", "length"]
        assert [row["read_index"] for row in json_output] == list(
            range(1, len(json_output) + 1)
        )

    def test_constant_sites_verbose_contract(self):
        target = sample_files / "simple.fa"

        tsv_output = _run_and_parse("constant_sites", target, "tsv", ["-v"])
        json_output = _run_and_parse("constant_sites", target, "json", ["-v"])
        yaml_output = _run_and_parse("constant_sites", target, "yaml", ["-v"])

        tsv_rows = []
        for row in tsv_output.splitlines():
            site_index, classification = row.split("\t")
            tsv_rows.append((int(site_index), classification))

        json_rows = [(row["site_index"], row["classification"]) for row in json_output]
        yaml_rows = [(row["site_index"], row["classification"]) for row in yaml_output]

        assert tsv_rows == json_rows
        assert tsv_rows == yaml_rows
        assert list(json_output[0].keys()) == ["site_index", "classification"]
        assert list(yaml_output[0].keys()) == ["site_index", "classification"]

    def test_rscu_contract(self, tmp_path):
        fasta = tmp_path / "cds.fa"
        translation_table = tmp_path / "tt.txt"
        fasta.write_text(">g1\nATGTTTTTC\n>g2\nATGTTTTTC\n")
        translation_table.write_text("AUG M\nUUU F\nUUC F\n")

        tsv_output = _run_and_parse(
            "relative_synonymous_codon_usage",
            fasta,
            "tsv",
            ["-tt", str(translation_table)],
        )
        json_output = _run_and_parse(
            "relative_synonymous_codon_usage",
            fasta,
            "json",
            ["-tt", str(translation_table)],
        )
        yaml_output = _run_and_parse(
            "relative_synonymous_codon_usage",
            fasta,
            "yaml",
            ["-tt", str(translation_table)],
        )

        tsv_map = {}
        for row in tsv_output.splitlines():
            codon, rscu = row.split("\t")
            tsv_map[codon] = float(rscu)

        json_map = {row["codon"]: row["rscu"] for row in json_output}
        yaml_map = {row["codon"]: row["rscu"] for row in yaml_output}

        assert tsv_map == json_map
        assert tsv_map == yaml_map
        assert list(json_output[0].keys()) == ["codon", "rscu"]
        assert list(yaml_output[0].keys()) == ["codon", "rscu"]

    def test_gene_wise_rscu_contract(self, tmp_path):
        fasta = tmp_path / "cds.fa"
        translation_table = tmp_path / "tt.txt"
        fasta.write_text(">g2\nATGTTTTTC\n>g1\nATGTTTTTC\n")
        translation_table.write_text("AUG M\nUUU F\nUUC F\n")

        tsv_output = _run_and_parse(
            "gene_wise_relative_synonymous_codon_usage",
            fasta,
            "tsv",
            ["-tt", str(translation_table)],
        )
        json_output = _run_and_parse(
            "gene_wise_relative_synonymous_codon_usage",
            fasta,
            "json",
            ["-tt", str(translation_table)],
        )
        yaml_output = _run_and_parse(
            "gene_wise_relative_synonymous_codon_usage",
            fasta,
            "yaml",
            ["-tt", str(translation_table)],
        )

        tsv_rows = {}
        for row in tsv_output.splitlines():
            fields = row.split("\t")
            tsv_rows[fields[0]] = {
                "mean_rscu": float(fields[1]),
                "median_rscu": float(fields[2]),
                "stddev_rscu": float(fields[3]),
            }

        json_rows = _rows_to_map(json_output, "gene_id")
        yaml_rows = _rows_to_map(yaml_output, "gene_id")

        assert tsv_rows == json_rows
        assert tsv_rows == yaml_rows
        assert list(json_output[0].keys()) == [
            "gene_id",
            "mean_rscu",
            "median_rscu",
            "stddev_rscu",
        ]
        assert list(yaml_output[0].keys()) == [
            "gene_id",
            "mean_rscu",
            "median_rscu",
            "stddev_rscu",
        ]
        assert [row["gene_id"] for row in json_output] == ["g1", "g2"]

    def test_fastq_read_lengths_summary_contract(self):
        target = sample_files / "DRR284700_1_subset.fastq"

        tsv_output = _run_and_parse("fastq_read_lengths", target, "tsv")
        json_output = _run_and_parse("fastq_read_lengths", target, "json")
        yaml_output = _run_and_parse("fastq_read_lengths", target, "yaml")

        mean, stdev = tsv_output.split(" +/- ")
        tsv_summary = {"mean": float(mean), "stdev": float(stdev)}

        assert tsv_summary == json_output
        assert tsv_summary == yaml_output
        assert list(json_output.keys()) == ["mean", "stdev"]
        assert list(yaml_output.keys()) == ["mean", "stdev"]

    def test_trim_se_fastq_summary_contract(self, tmp_path):
        fastq = tmp_path / "reads.fq"
        output_file = tmp_path / "reads.trimmed.fq"
        fastq.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGA\n+\nIIII\n")

        tsv_output = _run_and_parse(
            "trim_se_fastq",
            fastq,
            "tsv",
            ["-o", str(output_file), "-m", "0", "-l", "1"],
        )
        json_output = _run_and_parse(
            "trim_se_fastq",
            fastq,
            "json",
            ["-o", str(output_file), "-m", "0", "-l", "1"],
        )
        yaml_output = _run_and_parse(
            "trim_se_fastq",
            fastq,
            "yaml",
            ["-o", str(output_file), "-m", "0", "-l", "1"],
        )

        tsv_summary = {}
        for row in tsv_output.splitlines():
            label, value = row.split(": ")
            key = label.lower().replace(" ", "_")
            tsv_summary[key] = int(value)

        assert tsv_summary == json_output
        assert tsv_summary == yaml_output
        assert list(json_output.keys()) == [
            "reads_processed",
            "reads_kept",
            "reads_removed",
        ]
        assert list(yaml_output.keys()) == [
            "reads_processed",
            "reads_kept",
            "reads_removed",
        ]
