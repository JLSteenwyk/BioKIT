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
                line = line[2:]
                key, value = line.split(": ", 1)
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


def _run_cli(args: list[str]) -> str:
    with patch("builtins.print") as mocked_print:
        with patch.object(sys, "argv", args):
            Biokit()

    return "\n".join(str(call.args[0]) for call in mocked_print.call_args_list)


def _run_and_parse(command: str, target_file: Path, output_format: str, extra_args=None):
    extra_args = extra_args or []
    output = _run_cli(
        ["biokit", command, str(target_file), *extra_args, "--format", output_format]
    )

    if output_format == "json":
        return json.loads(output)
    if output_format == "yaml":
        return _parse_yaml_simple(output)
    return output


@pytest.mark.integration
class TestOutputFormatMatrix(object):
    def test_sequence_length_matrix(self):
        target = sample_files / "EOG091N44MS.fa"

        tsv_output = _run_and_parse("sequence_length", target, "tsv")
        json_output = _run_and_parse("sequence_length", target, "json")
        yaml_output = _run_and_parse("sequence_length", target, "yaml")

        tsv_rows = []
        for row in tsv_output.splitlines():
            identifier, length = row.split("\t")
            tsv_rows.append({"id": identifier, "length": int(length)})
        tsv_rows.sort(key=lambda row: row["id"])

        assert tsv_rows == json_output
        assert tsv_rows == yaml_output

    def test_genome_assembly_metrics_matrix(self):
        target = sample_files / "simple.fa"

        tsv_output = _run_and_parse("genome_assembly_metrics", target, "tsv")
        json_output = _run_and_parse("genome_assembly_metrics", target, "json")
        yaml_output = _run_and_parse("genome_assembly_metrics", target, "yaml")

        tsv_map = {}
        for row in tsv_output.splitlines():
            value, label = row.split("\t")
            tsv_map[GENOME_ASSEMBLY_LABEL_TO_KEY[label]] = _parse_scalar(value)

        assert tsv_map == json_output
        assert tsv_map == yaml_output

    def test_fastq_read_lengths_verbose_matrix(self):
        target = sample_files / "DRR284700_1_subset.fastq"

        tsv_output = _run_and_parse("fastq_read_lengths", target, "tsv", ["-v"])
        json_output = _run_and_parse("fastq_read_lengths", target, "json", ["-v"])
        yaml_output = _run_and_parse("fastq_read_lengths", target, "yaml", ["-v"])

        tsv_lengths = [int(line) for line in tsv_output.splitlines()]
        json_lengths = [row["length"] for row in json_output]
        yaml_lengths = [row["length"] for row in yaml_output]

        assert tsv_lengths == json_lengths
        assert tsv_lengths == yaml_lengths
