#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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

Count, normalize and embed k-mers given nucleotide sequences
"""

import pytest
import os

import pandas
from unittest.mock import patch, MagicMock

from io import StringIO, BytesIO
from autometa.common import kmers
from autometa.common.exceptions import KmerFormatError, KmerEmbeddingError

assembly = os.path.join("tests", "data", "metagenome.fna")


def test_revcomp():
    string = "AAATTGCGCCCCCG"
    rev_string = kmers._revcomp(string)
    assert rev_string == "CGGGGGCGCAATTT"
    assert type(rev_string) is str
    assert len(rev_string) == 14


def test_init_kmers():
    unique_kmers_5 = kmers.init_kmers(kmer_size=5)
    assert type(unique_kmers_5) is dict
    assert len(unique_kmers_5) == 512
    unique_kmers_3 = kmers.init_kmers(kmer_size=3.0)
    assert type(unique_kmers_3) is dict
    assert len(unique_kmers_3) == 32
    unique_kmers_5 = kmers.init_kmers(kmer_size="5")
    assert type(unique_kmers_5) is dict
    assert len(unique_kmers_5) == 512
    unique_kmers_3 = kmers.init_kmers(kmer_size="3.0")
    assert type(unique_kmers_3) is dict
    assert len(unique_kmers_3) == 32
    with pytest.raises(TypeError):
        kmers.init_kmers(kmer_size=5.5)
    with pytest.raises(TypeError):
        kmers.init_kmers(kmer_size="5.5")


@patch("os.path.getsize", return_value=2 * 1024 * 1024)
def test_load(patched_file_size):
    kmer_table_fpath = MagicMock(
        return_value="path_to_kmer_frequency_table",
        name="path to input kmer frequency table",
    )
    with pytest.raises(FileNotFoundError):
        kmers.load("kmer_fpath")
    with pytest.raises(KmerFormatError):
        kmers.load(kmer_table_fpath)
        assert patched_file_size.called is True


# def test_seq_counter():
#     ref_kmers = kmers.init_kmers()
#     kmer_counts = kmers.seq_counter(
#         assembly=assembly, ref_kmers=ref_kmers, verbose=True
#     )
#     assert type(kmer_counts) is dict


@patch("test_kmers.kmers.normalize")
@patch("test_kmers.kmers.seq_counter")
@patch("test_kmers.kmers.mp_counter")
def test_count(pacthed_map_counter, pacthed_seq_counter, patched_normalize):
    with pytest.raises(TypeError):
        kmers.count(assembly=assembly, kmer_size=5.5)
    with pytest.raises(TypeError):
        kmers.count(assembly=assembly, kmer_size="5.5")

    out_df = kmers.count(assembly=assembly, kmer_size=5.0)
    assert type(out_df) is pd.DataFrame
    assert out_df.index.name == "contig"
    assert len(list(out_df)) == 512

    out_df = kmers.count(assembly=assembly, kmer_size=5, multiprocess=False)
    assert type(out_df) is pd.DataFrame
    assert out_df.index.name == "contig"
    assert len(list(out_df)) == 512

    out_df = kmers.count(assembly=assembly, kmer_size="5.0", normalized=True)

    assert pacthed_map_counter.call_count == 2
    assert pacthed_seq_counter.call_count == 1
    assert patched_normalize.call_count == 1


# @patch("pd.DataFrame.empty", return_value=False)
# def test_embed():
#     kmer_fpath = MagicMock(
#         return_value="path_to_kmers_file", name="path to kmer frequency table",
#     )

#     with pytest.raises(KmerEmbeddingError):
#         kmers.embed(kmers=None, embedded=None)

#     # with pytest.raises(ValueError):
#     #     kmers.embed(kmers=kmer_fpath)
#     # with patch(
#     #     "pd.read_csv",
#     #     return_value={"contig1": ["count1", "count2"], "contig2": ["count1", "count2"]},
#     # ):
#     #     kmers.embed(kmers=kmer_fpath)
#     test_df = pd.DataFrame(
#         {
#             "contig": [
#                 "NODE_1_length_1389215_cov_225.275",
#                 "NODE_2_length_1166739_cov_224.155",
#                 "NODE_3_length_1063064_cov_225.095",
#                 "NODE_4_length_1031470_cov_223.812",
#                 "NODE_5_length_937195_cov_225.122",
#             ],
#             "AAAAA": [170, 1688, 143, 347, 156],
#             "AAAAG": [412, 1354, 345, 331, 317],
#             "AAATG": [536, 1513, 286, 325, 378],
#         }
#     )
#     output = BytesIO()
#     writer = pd.ExcelWriter(output, engine="xlsxwriter")
#     test_df.to_excel(writer, sheet_name="Sheet1", index=False)
#     writer.save()
#     output.seek(0)

#     kmer_fpath = MagicMock(return_value=output, name="path to kmer frequency table",)
#     kmers.embed(kmers=test_df)


# @patch("os.path.getsize", return_value=2 * 1024 * 1024)
# @patch("pandas.DataFrame.empty", return_value=False)
# def test_embed(patched_df_empty, patched_file_size):
#     with pytest.raises(KmerEmbeddingError):
#         kmers.embed(kmers=None, embedded=None)

#     kmer_fpath = MagicMock(
#         return_value="path_to_kmers_file", name="path to kmer frequency table"
#     )

#     kmers.embed(kmers=kmer_fpath)
