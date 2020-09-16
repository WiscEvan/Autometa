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

import pandas as pd
from unittest.mock import patch, MagicMock

from autometa.common import kmers
from autometa.common.exceptions import TableFormatError

assembly = os.path.join("tests", "data", "metagenome.fna")
test_df = pd.DataFrame(
    {
        "contig": [
            "NODE_1_length_1389215_cov_225.275",
            "NODE_2_length_1166739_cov_224.155",
            "NODE_3_length_1063064_cov_225.095",
            "NODE_4_length_1031470_cov_223.812",
            "NODE_5_length_937195_cov_225.122",
        ],
        "AAAAA": [170, 1688, 143, 347, 156],
        "AAAAG": [412, 1354, 345, 331, 317],
        "AAATG": [536, 1513, 286, 325, 378],
    }
)
test_df.set_index("contig", inplace=True)


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
    with pytest.raises(TableFormatError):
        kmers.load(kmer_table_fpath)
    assert patched_file_size.called is True


@patch("autometa.common.kmers.seq_counter")
@patch("autometa.common.kmers.mp_counter")
def test_count(pacthed_mp_counter, pacthed_seq_counter):
    with pytest.raises(TypeError):
        kmers.count(assembly=assembly, size=5.5)
    with pytest.raises(TypeError):
        kmers.count(assembly=assembly, size="5.5")

    out_df = kmers.count(assembly=assembly, size=5.0)
    assert type(out_df) is pd.DataFrame
    assert out_df.index.name == "contig"
    assert out_df.shape == (0, 512)

    out_df = kmers.count(assembly=assembly, size=5, multiprocess=False)
    assert type(out_df) is pd.DataFrame
    assert out_df.index.name == "contig"
    assert out_df.shape == (0, 512)

    assert pacthed_mp_counter.call_count == 1
    assert pacthed_seq_counter.call_count == 1


@patch("autometa.common.kmers.autometa_clr", side_effect=kmers.autometa_clr)
def test_normalize(pacthed_autometa_clr, tmpdir):
    out_fpath = MagicMock(
        return_value="path_to_output_table", name="path to input kmer normalised table",
    )
    norm_df = kmers.normalize(df=test_df)
    assert type(norm_df) is pd.DataFrame
    assert norm_df.index.name == "contig"
    assert list(norm_df) == ["AAAAA", "AAAAG", "AAATG"]
    assert round((norm_df.iloc[2, 2]), 2) == 0.17
    assert round((norm_df.iloc[3, 1]), 2) == -0.01
    assert norm_df.shape == (5, 3)

    norm_df = kmers.normalize(df=test_df, out=tmpdir.join("out_table"))
    assert type(norm_df) is pd.DataFrame
    assert norm_df.index.name == "contig"
    assert list(norm_df) == ["AAAAA", "AAAAG", "AAATG"]
    assert round((norm_df.iloc[2, 2]), 2) == 0.17
    assert round((norm_df.iloc[3, 1]), 2) == -0.01
    assert norm_df.shape == (5, 3)

    norm_df = kmers.normalize(df=test_df, out=out_fpath, force=True)
    assert type(norm_df) is pd.DataFrame
    assert norm_df.index.name == "contig"
    assert list(norm_df) == ["AAAAA", "AAAAG", "AAATG"]
    assert round((norm_df.iloc[2, 2]), 2) == 0.17
    assert round((norm_df.iloc[3, 1]), 2) == -0.01
    assert norm_df.shape == (5, 3)

    norm_df = kmers.normalize(df=test_df, method="ilr")
    assert type(norm_df) is pd.DataFrame
    assert norm_df.index.name == "contig"
    assert list(norm_df) == [0, 1]
    assert round((norm_df.iloc[2, 1]), 2) == -0.21
    assert round((norm_df.iloc[3, 0]), 2) == 0.03
    assert norm_df.shape == (5, 2)

    norm_df = kmers.normalize(df=test_df, method="clr")
    assert type(norm_df) is pd.DataFrame
    assert norm_df.index.name == "contig"
    assert list(norm_df) == [0, 1, 2]
    assert round((norm_df.iloc[2, 2]), 2) == 0.17
    assert round((norm_df.iloc[3, 1]), 2) == -0.01
    assert norm_df.shape == (5, 3)

    with pytest.raises(ValueError):
        kmers.normalize(df=test_df, method="invalid_method")

    assert pacthed_autometa_clr.call_count == 3


@patch("os.path.getsize", return_value=2 * 1024 * 1024)
def test_embed(patched_file_size):
    # spec is needed as it makes the instance of MagicMock to str, but if we use spec os.path.exists gives us False, thus we need to patch that
    #  See https://stackoverflow.com/a/11283173/12671809
    kmer_fpath = MagicMock(
        return_value="path_to_kmers_file",
        name="path to kmer frequency table",
        spec="path_to_kmers_file",
    )
    out_fpath = MagicMock(
        return_value="path_to_output_file", name="path to output table",
    )

    with patch("os.path.exists", return_value=True) as mock_path:
        with pytest.raises(TableFormatError):
            kmers.embed(kmers=kmer_fpath)
    assert mock_path.called is True

    with pytest.raises(TypeError):
        kmers.embed(kmers="invalid_fpath")

    normalized_test_df = kmers.normalize(df=test_df)
    with pytest.raises(TableFormatError):
        kmers.embed(kmers=normalized_test_df, out=out_fpath)

    empty_df = pd.DataFrame({})
    with pytest.raises(FileNotFoundError):
        kmers.embed(kmers=empty_df, out=out_fpath, force=True)

    with pytest.raises(ValueError):
        kmers.embed(kmers=normalized_test_df, method="invalid_method")
    # Need to test if it raises value error while doing dimension reduction but not sure how to fail this step
    # with pytest.raises(ValueError):
    #     kmers.embed(kmers=normalized_test_df, method="bhsne", embed_dimensions=10)

    kmers.embed(kmers=normalized_test_df, method="sksne", embed_dimensions=10)

    embed_df = kmers.embed(kmers=normalized_test_df, method="bhsne")
    assert type(embed_df) is pd.DataFrame
    assert embed_df.index.name == "contig"
    assert list(embed_df) == ["x", "y"]
    assert round((embed_df.iloc[2, 1]), 2) == -618.96
    assert round((embed_df.iloc[4, 1]), 2) == 3395.94
    assert embed_df.shape == (5, 2)

    embed_df = kmers.embed(kmers=normalized_test_df, out=out_fpath, force=True)
    assert type(embed_df) is pd.DataFrame
    assert embed_df.index.name == "contig"
    assert list(embed_df) == ["x", "y"]
    assert round((embed_df.iloc[2, 1]), 2) == -618.96
    assert round((embed_df.iloc[4, 1]), 2) == 3395.94
    assert embed_df.shape == (5, 2)

    embed_df = kmers.embed(kmers=normalized_test_df, embed_dimensions=2, method="sksne")
    # Coversion to float64 is needed for sksne and umap as in these cases the DataFrame returns values which are in float32, and thus can't be compared with normal float.
    embed_df = embed_df.astype("float64")
    assert type(embed_df) is pd.DataFrame
    assert embed_df.index.name == "contig"
    assert list(embed_df) == ["x", "y"]
    assert round((embed_df.iloc[2, 1]), 2) == -667.56
    assert round((embed_df.iloc[4, 1]), 2) == -957.08
    assert embed_df.shape == (5, 2)

    embed_df_10_dim = kmers.embed(
        kmers=normalized_test_df, embed_dimensions=10, method="sksne"
    )
    embed_df_10_dim = embed_df_10_dim.astype("float64")
    assert type(embed_df_10_dim) is pd.DataFrame
    assert embed_df_10_dim.index.name == "contig"
    assert list(embed_df_10_dim) == ["x", "y", "z"]
    assert round((embed_df_10_dim.iloc[2, 1]), 2) == 412.51
    assert round((embed_df_10_dim.iloc[4, 1]), 2) == -911.67
    assert embed_df_10_dim.shape == (5, 3)

    embed_df = kmers.embed(kmers=normalized_test_df, embed_dimensions=3, method="sksne")
    embed_df = embed_df.astype("float64")
    assert type(embed_df) is pd.DataFrame
    assert embed_df.index.name == "contig"
    assert list(embed_df) == ["x", "y", "z"]
    assert round((embed_df.iloc[2, 1]), 2) == 412.51
    assert round((embed_df.iloc[4, 1]), 2) == -911.67
    pd.testing.assert_frame_equal(embed_df_10_dim, embed_df)
    assert embed_df.shape == (5, 3)

    embed_df = kmers.embed(kmers=normalized_test_df, embed_dimensions=2, method="umap")
    embed_df = embed_df.astype("float64")
    assert type(embed_df) is pd.DataFrame
    assert embed_df.index.name == "contig"
    assert list(embed_df) == ["x", "y"]
    assert round((embed_df.iloc[2, 1]), 2) == -5.2
    assert round((embed_df.iloc[4, 1]), 2) == -4.28
    assert embed_df.shape == (5, 2)

    assert patched_file_size.call_count == 4
