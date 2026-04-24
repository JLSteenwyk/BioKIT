from argparse import Namespace

import pytest

from biokit.services.text.sample_sequences import SampleSequences


def _args(fasta, number=None, percent=None, seed=None, output=None):
    return Namespace(
        fasta=fasta,
        number=number,
        percent=percent,
        seed=seed,
        output=output,
    )


def _write_fasta(tmp_path, n=10):
    fa = tmp_path / "seqs.fa"
    with open(fa, "w") as fh:
        for i in range(n):
            fh.write(f">seq{i} some description {i}\nATCG{i:04d}\n")
    return str(fa)


def _parse_output_headers(out):
    return [line[1:].split()[0] for line in out.splitlines() if line.startswith(">")]


def test_sample_sequences_number(tmp_path, capsys):
    fa = _write_fasta(tmp_path, n=10)
    SampleSequences(_args(fa, number=3, seed=42)).run()
    headers = _parse_output_headers(capsys.readouterr().out)
    assert len(headers) == 3
    assert all(h.startswith("seq") for h in headers)
    assert len(set(headers)) == 3  # no duplicates


def test_sample_sequences_percent(tmp_path, capsys):
    fa = _write_fasta(tmp_path, n=20)
    SampleSequences(_args(fa, percent=25, seed=0)).run()
    headers = _parse_output_headers(capsys.readouterr().out)
    assert len(headers) == 5  # 25% of 20


def test_sample_sequences_default_percent(tmp_path, capsys):
    fa = _write_fasta(tmp_path, n=50)
    SampleSequences(_args(fa, seed=0)).run()
    headers = _parse_output_headers(capsys.readouterr().out)
    # default percent = 10% of 50
    assert len(headers) == 5


def test_sample_sequences_seed_reproducible(tmp_path, capsys):
    fa = _write_fasta(tmp_path, n=50)

    SampleSequences(_args(fa, number=10, seed=123)).run()
    first = capsys.readouterr().out

    SampleSequences(_args(fa, number=10, seed=123)).run()
    second = capsys.readouterr().out

    assert first == second


def test_sample_sequences_different_seeds_differ(tmp_path, capsys):
    fa = _write_fasta(tmp_path, n=50)

    SampleSequences(_args(fa, number=10, seed=1)).run()
    first = capsys.readouterr().out

    SampleSequences(_args(fa, number=10, seed=2)).run()
    second = capsys.readouterr().out

    assert first != second


def test_sample_sequences_preserves_full_header(tmp_path, capsys):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">seq1 hello world\nATCG\n>seq2 foo bar\nGGGG\n")
    SampleSequences(_args(str(fa), number=2, seed=0)).run()
    out = capsys.readouterr().out
    assert ">seq1 hello world" in out
    assert ">seq2 foo bar" in out


def test_sample_sequences_output_file(tmp_path):
    fa = _write_fasta(tmp_path, n=10)
    out_path = tmp_path / "sample.fa"
    SampleSequences(_args(fa, number=3, seed=0, output=str(out_path))).run()
    content = out_path.read_text()
    assert content.count(">") == 3


def test_sample_sequences_number_equals_total(tmp_path, capsys):
    fa = _write_fasta(tmp_path, n=5)
    SampleSequences(_args(fa, number=5, seed=0)).run()
    headers = _parse_output_headers(capsys.readouterr().out)
    assert len(headers) == 5


def test_sample_sequences_zero_number(tmp_path, capsys):
    fa = _write_fasta(tmp_path, n=10)
    SampleSequences(_args(fa, number=0, seed=0)).run()
    out = capsys.readouterr().out
    assert out == ""


def test_sample_sequences_rejects_both_number_and_percent(tmp_path):
    fa = _write_fasta(tmp_path, n=10)
    with pytest.raises(ValueError):
        SampleSequences(_args(fa, number=3, percent=10))


def test_sample_sequences_rejects_negative_number(tmp_path):
    fa = _write_fasta(tmp_path, n=10)
    with pytest.raises(ValueError):
        SampleSequences(_args(fa, number=-1))


def test_sample_sequences_rejects_percent_out_of_range(tmp_path):
    fa = _write_fasta(tmp_path, n=10)
    with pytest.raises(ValueError):
        SampleSequences(_args(fa, percent=150))


def test_sample_sequences_rejects_number_greater_than_total(tmp_path):
    fa = _write_fasta(tmp_path, n=5)
    with pytest.raises(ValueError):
        SampleSequences(_args(fa, number=10, seed=0)).run()
