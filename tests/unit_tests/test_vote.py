#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

Script to test autometa/taxonomy/vote.py
"""

import argparse
import pytest

from autometa.taxonomy import vote

import pandas as pd
from autometa.common.exceptions import TableFormatError
from Bio import SeqIO


@pytest.fixture(name="blastp")
def fixture_blastp(variables, tmp_path):
    """diamond blastp output fixture

    Parameters
    ----------
    variables : dict
        Imported from test_data.json (See pytest.ini for reference to file path).
        test_data.json is generated by make_test_data.py
    tmp_path : Pathlib.Path
        temporary path constructed by pytest

    Returns
    -------
    str
        diamond blastp output file path written from `variables["taxonomy"]`
    """
    vote_test_data = variables["taxonomy"]
    blastp_fpath = tmp_path / "blastp.tsv"
    df = pd.read_json(vote_test_data["blastp"])
    df.to_csv(blastp_fpath, sep="\t", index=False, header=False)
    return str(blastp_fpath)


@pytest.fixture(name="prot_orfs")
def fixture_prot_orfs(variables, tmp_path):
    """Amino-acid orfs corresponding to diamond blastp output table.
    This is required for translation of query ORFs back to contigs.

    Parameters
    ----------
    variables : dict
        Imported from test_data.json (See pytest.ini for reference to file path).
        test_data.json is generated by make_test_data.py
    tmp_path : Pathlib.Path
        temporary path constructed by pytest

    Returns
    -------
    str
        amino-acid ORFs file path
    """
    vote_test_data = variables["taxonomy"]
    records = vote_test_data["prot_orfs"]
    lines = ""
    for record, seq in records.items():
        lines += f"{record}\n{seq}\n"
    fpath = tmp_path / "orfs.faa"
    with open(fpath, "w") as fh:
        fh.write(lines)
    return str(fpath)


@pytest.fixture(name="votes", scope="module")
def fixture_votes():
    votes = [
        {"contig": "NODE_1_length_1389215_cov_225.275", "taxid": 373},
        {"contig": "NODE_2_length_1166739_cov_224.155", "taxid": 60890},
    ]
    return pd.DataFrame(votes).set_index("contig")


@pytest.fixture(name="votes_fpath", scope="module")
def fixture_votes_fpath(votes, tmp_path_factory):
    fpath = tmp_path_factory.mktemp("vote") / "votes.tsv"
    votes.to_csv(fpath, sep="\t", index=True, header=True)
    return fpath


def test_add_ranks(ncbi, votes, tmp_path):
    out = tmp_path / "taxonomy.ranks_added.tsv"
    df = vote.add_ranks(df=votes, ncbi=ncbi)
    assert df.shape == (2, 8)
    assert df.index.name == "contig"
    canonical_ranks = {rank for rank in ncbi.CANONICAL_RANKS if rank != "root"}
    for canonical_rank in canonical_ranks:
        assert canonical_rank in df.columns


@pytest.mark.slow
@pytest.mark.skip
def test_vote_assign(blastp, ncbi_dir, prot_orfs, tmp_path):
    out = tmp_path / "votes.tsv"
    votes = vote.assign(
        out=out,
        prot_orfs=prot_orfs,
        blast=blastp,
        ncbi_dir=ncbi_dir,
    )
    assert isinstance(votes, pd.DataFrame)
    assert votes.index.name == "contig"
    assert "taxid" in votes.columns


def test_get(ncbi, votes_fpath):
    df = vote.get(
        filepath_or_dataframe=votes_fpath,
        kingdom="bacteria",
        ncbi=ncbi,
    )
    # canonical ranks should have been added to table if they were not already in place.
    assert df.shape == (2, 8)


def test_get_none_recovered(ncbi, votes_fpath):
    with pytest.raises(KeyError):
        vote.get(
            filepath_or_dataframe=votes_fpath,
            kingdom="archaea",
            ncbi=ncbi,
        )


def test_get_empty_votes(ncbi_dir, tmp_path):
    fpath = tmp_path / "votes.tsv"
    with pytest.raises(FileNotFoundError):
        vote.get(
            filepath_or_dataframe=fpath,
            kingdom="archaea",
            ncbi=ncbi_dir,
        )


def test_get_superkingdom_not_in_columns(monkeypatch, ncbi, votes, tmp_path):
    def return_df(*args, **kwargs):
        return votes

    fpath = tmp_path / "votes.tsv"
    votes.to_csv(fpath, sep="\t", index=True, header=True)
    monkeypatch.setattr(vote, "add_ranks", return_df, raising=True)
    with pytest.raises(TableFormatError):
        vote.get(
            filepath_or_dataframe=fpath,
            kingdom="archaea",
            ncbi=ncbi,
        )


@pytest.fixture(name="ranks_added_votes", scope="module")
def fixture_ranks_added_votes(votes_fpath, ncbi):
    return vote.get(
        filepath_or_dataframe=votes_fpath,
        kingdom="bacteria",
        ncbi=ncbi,
    )


@pytest.mark.parametrize(
    "rank,num_expected", [("superkingdom", 1), ("noncanonical_rank", 0), ("order", 2)]
)
def test_write_ranks(monkeypatch, tmp_path, ranks_added_votes, rank, num_expected):
    dirpath = tmp_path / "metabins"
    dirpath.mkdir()
    assembly = dirpath / "assembly.fna"
    assembly.write_text("records")
    if rank == "noncanonical_rank":
        with pytest.raises(ValueError):
            vote.write_ranks(
                taxonomy=ranks_added_votes, assembly=assembly, outdir=dirpath, rank=rank
            )
        return

    class MockedRecord:
        def __init__(self, id):
            self.id = id

    def return_mock_assembly_records(*args, **kwargs):
        return [
            MockedRecord(id="NODE_1_length_1389215_cov_225.275"),
            MockedRecord(id="NODE_2_length_1166739_cov_224.155"),
        ]

    def return_mock_write(*args, **kwargs):
        return 2

    monkeypatch.setattr(SeqIO, "parse", return_mock_assembly_records, raising=True)
    monkeypatch.setattr(SeqIO, "write", return_mock_write)
    fpaths = vote.write_ranks(
        taxonomy=ranks_added_votes, assembly=assembly, outdir=dirpath, rank=rank
    )
    assert len(fpaths) == num_expected


def test_write_ranks_no_assembly(tmp_path, ranks_added_votes):
    dirpath = tmp_path / "metabins"
    dirpath.mkdir()
    assembly = dirpath / "assembly.fna"
    with pytest.raises(FileNotFoundError):
        vote.write_ranks(
            taxonomy=ranks_added_votes,
            assembly=assembly,
            outdir=dirpath,
            rank="superkingdom",
        )


def test_write_ranks_no_taxonomy_columns(tmp_path, votes):
    dirpath = tmp_path / "metabins"
    dirpath.mkdir()
    assembly = dirpath / "assembly.fna"
    with pytest.raises(KeyError):
        vote.write_ranks(
            taxonomy=votes,
            assembly=assembly,
            outdir=dirpath,
            rank="superkingdom",
        )


@pytest.mark.slow
@pytest.mark.skip
@pytest.mark.entrypoint
def test_vote_main(monkeypatch, ncbi_dir, tmp_path):
    outdir = tmp_path / "outdir"
    outdir.mkdir()
    taxonomy = outdir / "taxonomy.tsv"
    assembly = outdir / "assembly.fna"
    assembly.write_text("records")
    lca_fpath = outdir / "lca.tsv"
    lca_fpath = str(lca_fpath)

    class MockArgs:
        def __init__(self):
            self.input = taxonomy
            self.output = outdir
            self.assembly = assembly
            self.split_rank_and_write = "superkingdom"
            self.ncbi = ncbi_dir

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    def return_mock_parser(*args, **kwargs):
        return MockParser()

    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)
    vote.main()
    assert taxonomy.exists()
