import pytest
import os
import argparse

import pandas as pd
from Bio import SeqIO

from autometa.common import kmers
from unittest.mock import patch, MagicMock
from autometa.common.exceptions import TableFormatError


@pytest.fixture(name="assembly", scope="module")
def fixture_assembly(variables, tmp_path_factory):
    kmer_test_data = variables["metagenome"]
    records = kmer_test_data["assembly"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    outdir = tmp_path_factory.mktemp("kmers")
    fpath = outdir / "metagenome.fna"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="counts", scope="module")
def fixture_counts(variables):
    kmer_test_data = variables["kmers"]
    df = pd.read_json(kmer_test_data["counts"])
    # kmer size is 5 (b/c this is the default).
    df.set_index("contig", inplace=True)
    return df


@pytest.fixture(name="counts_fpath", scope="module")
def fixture_counts_fpath(counts, tmp_path_factory):
    fpath = tmp_path_factory.mktemp("kmers") / "counts.tsv"
    counts.to_csv(fpath, sep="\t", index=True, header=True)
    return fpath.as_posix()


@pytest.fixture(name="norm_df", scope="module")
def fixture_norm_df(variables):
    kmer_test_data = variables["kmers"]
    df = pd.read_json(kmer_test_data["am_clr_normalized_counts"])
    df.set_index("contig", inplace=True)
    return df


@pytest.fixture(name="invalid_df_fpath")
def fixture_df_without_contig_index_(tmp_path):
    invalid_dict = {
        "column1": ["invalid_contig_1", "invalid_contig_2", "invalid_contig_3"],
        "column2": ["invalid_marker1", "invalid_marker2", "invalid_marker3"],
    }
    df = pd.DataFrame(invalid_dict)
    df_fpath = tmp_path / "invalid_df.tsv"
    df.to_csv(df_fpath)
    return df_fpath.as_posix()


def test_kmer_load(counts_fpath):
    df = kmers.load(kmers_fpath=counts_fpath)
    assert not df.empty
    assert df.index.name == "contig"


def test_kmer_load_FileNotFoundError():
    with pytest.raises(FileNotFoundError):
        kmers.load(kmers_fpath="Invalid_fpath")


def test_kmer_load_TableFormatError(invalid_df_fpath):
    with pytest.raises(TableFormatError):
        kmers.load(invalid_df_fpath)


@pytest.mark.parametrize("multiprocess", [True, False])
def test_count(assembly, multiprocess, tmp_path):
    out = tmp_path / "kmers.tsv"
    size = 5
    force = False
    df = kmers.count(
        assembly=assembly, size=size, out=out, force=force, multiprocess=multiprocess,
    )
    assert df.shape[1] == 4 ** size / 2
    assert df.index.name == "contig"
    assert out.exists()


@pytest.mark.parametrize("force", [True, False])
def test_count_out_exists(assembly, counts, force, tmp_path):
    out = tmp_path / "kmers.tsv"
    counts.to_csv(out, sep="\t", index=True, header=True)
    size = 5
    df = kmers.count(
        assembly=assembly, size=size, out=out, force=force, multiprocess=True,
    )
    assert df.shape[1] == 4 ** size / 2
    assert df.index.name == "contig"
    assert out.exists()


def test_count_wrong_size(assembly, tmp_path):
    size = 5.5
    with pytest.raises(TypeError):
        kmers.count(assembly=assembly, size=size)


@pytest.mark.parametrize("method", ["am_clr", "clr", "ilr"])
def test_normalize(counts, method, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    force = False
    df = kmers.normalize(df=counts, method=method, out=out, force=force)
    if method in {"am_clr", "clr"}:
        assert df.shape == counts.shape
    else:
        # ILR will reduce the columns by one.
        assert df.shape[1] < counts.shape[1]
    assert out.exists()


@pytest.mark.parametrize("force", [True, False])
def test_normalize_out_exists(counts, norm_df, force, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    norm_df.to_csv(out, sep="\t", index=True, header=True)
    df = kmers.normalize(df=counts, method="am_clr", out=out, force=force)
    assert df.shape == counts.shape
    assert df.index.name == "contig"


def test_normalize_wrong_method(counts, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    with pytest.raises(ValueError):
        kmers.normalize(df=counts, method="am_ilr", out=out, force=False)


@pytest.mark.parametrize("method", ["bhsne", "sksne", "umap"])
def test_embed(norm_df, method, tmp_path):
    seed = 42
    out = tmp_path / "kmers.embed.tsv"
    force = False
    embed_dimensions = 2
    do_pca = True
    pca_dimensions = 3
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=force,
        embed_dimensions=embed_dimensions,
        do_pca=do_pca,
        pca_dimensions=pca_dimensions,
        method=method,
        seed=seed,
    )
    assert df.shape[1] == embed_dimensions


def test_embed_out_exists(norm_df, tmp_path):
    seed = 42
    out = tmp_path / "kmers.embed.tsv"
    force = False
    method = "bhsne"
    embed_dimensions = 2
    do_pca = True
    pca_dimensions = 3
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=force,
        embed_dimensions=embed_dimensions,
        do_pca=do_pca,
        pca_dimensions=pca_dimensions,
        method=method,
        seed=seed,
    )
    assert df.shape[1] == embed_dimensions
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=force,
        embed_dimensions=embed_dimensions,
        do_pca=do_pca,
        pca_dimensions=pca_dimensions,
        method=method,
        seed=seed,
    )


def test_embed_TableFormatError(invalid_df_fpath):
    with pytest.raises(TableFormatError):
        kmers.embed(kmers=invalid_df_fpath)


def test_embed_TypeError(tmp_path):
    kmer_fpath = tmp_path / "kmers.embed.tsv"
    with pytest.raises(TypeError):
        kmers.embed(kmers=kmer_fpath)


@patch("os.path.getsize", return_value=2 * 1024 * 1024)
def test_embed_FileNotFoundError(pacthed_file_size, tmp_path):
    empty_df = pd.DataFrame({})
    out = tmp_path / "kmers.embed.tsv"
    with pytest.raises(FileNotFoundError):
        kmers.embed(kmers=empty_df, out=out, force=True)


@pytest.fixture(name="mock_parser")
def fixture_mock_parser(
    assembly, counts, norm_df, monkeypatch, tmp_path,
):
    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockParseArgs:
        def __init__(self, assembly, counts, out, norm_df):
            self.fasta = assembly
            self.size = 5
            self.kmers = out
            self.force = True
            self.multiprocess = True
            self.cpus = 2
            self.normalized = True
            self.norm_method = "am_clr"
            self.normalized = out
            self.embedded = True
            self.embedded = out
            self.embed_method = "bhsne"
            self.embed_dimensions = 2
            self.do_pca = True
            self.pca_dimensions = 2
            self.seed = 42

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            out = tmp_path / "binning.tsv"
            return MockParseArgs(assembly, counts, out, norm_df,)

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)


def test_coverage_main(monkeypatch, mock_parser):
    kmers.main()
