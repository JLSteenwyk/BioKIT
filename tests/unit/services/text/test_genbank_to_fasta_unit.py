from argparse import Namespace
from pathlib import Path

from biokit.services.text.genbank_to_fasta import GenbankToFasta

SAMPLE = (
    Path(__file__).parent.parent.parent.parent
    / "sample_files"
    / "test_genbank.gb"
)


def _args(**overrides):
    defaults = dict(
        genbank=str(SAMPLE),
        feature_type=None,
        translate=False,
        output=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def test_genbank_to_fasta_whole_record(capsys):
    GenbankToFasta(_args()).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    assert lines[0].startswith(">TEST001.1")
    assert len(lines) == 2
    assert len(lines[1]) == 90


def test_genbank_to_fasta_cds_filter(capsys):
    GenbankToFasta(_args(feature_type="CDS")).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    # Two CDS features, 2 lines each
    assert len(lines) == 4
    assert ">TEST001_01" in lines[0]
    assert "CDS TEST001.1:1-30(+)" in lines[0]
    assert lines[1] == "ATGAAAATTTTATAGGCGCATTGAGCGCCT"
    assert ">TEST001_02" in lines[2]
    assert lines[3] == "ATGGCTCGCAACTAACGATGGGCGTTAAGG"


def test_genbank_to_fasta_rrna_filter(capsys):
    GenbankToFasta(_args(feature_type="rRNA")).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    assert len(lines) == 2
    assert ">TEST001_rrna_01" in lines[0]
    assert "rRNA" in lines[0]
    assert lines[1] == "CCATAAAATTGGCATTTTAAGCCTGCACAC"


def test_genbank_to_fasta_multiple_feature_types(capsys):
    GenbankToFasta(_args(feature_type="CDS,rRNA")).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    # 2 CDS + 1 rRNA = 3 records = 6 lines
    assert len(lines) == 6
    headers = [ln for ln in lines if ln.startswith(">")]
    assert len(headers) == 3


def test_genbank_to_fasta_translate_uses_qualifier(capsys):
    GenbankToFasta(_args(feature_type="CDS", translate=True)).run()
    out = capsys.readouterr().out
    lines = out.strip().split("\n")
    # Translation qualifiers from the sample file
    assert lines[1] == "MKIL*"
    assert lines[3] == "MARN*"


def test_genbank_to_fasta_output_file(tmp_path):
    out_path = tmp_path / "out.fa"
    GenbankToFasta(_args(output=str(out_path))).run()
    assert out_path.exists()
    content = out_path.read_text()
    assert content.startswith(">TEST001.1")
    assert "ATGAAAATTTTATAGG" in content


def test_genbank_to_fasta_case_insensitive_feature_type(capsys):
    GenbankToFasta(_args(feature_type="cds")).run()
    out = capsys.readouterr().out
    assert ">TEST001_01" in out
    assert ">TEST001_02" in out


def test_genbank_to_fasta_unknown_feature_type_yields_empty(capsys):
    GenbankToFasta(_args(feature_type="tRNA")).run()
    out = capsys.readouterr().out
    assert out == ""
