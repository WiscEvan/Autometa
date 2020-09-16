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

import os

from unittest.mock import mock_open, patch, MagicMock

from autometa.common.external import prodigal

assembly = os.path.join("tests", "data", "metagenome.fna")


contigs = [
    "NODE_11_length_590705_cov_223.126",
    "NODE_12_length_584723_cov_225.634",
    "NODE_13_length_533917_cov_222.896",
    "NODE_100_length_142487_cov_223.46",
    "NODE_110_length_132131_cov_223.987",
    "NODE_116_length_123186_cov_224.013",
    "NODE_119_length_119677_cov_223.947",
    "NODE_120_length_119298_cov_224.179",
    "NODE_132_length_112038_cov_223.935",
    "NODE_133_length_111688_cov_223.412",
    "NODE_136_length_111142_cov_224.73",
    "NODE_1074_length_13345_cov_224.047",
    "NODE_1076_length_13288_cov_223.745",
    "NODE_1080_length_13246_cov_222.788",
    "NODE_1086_length_13138_cov_224.524",
    "NODE_22_length_427389_cov_223.136",
    "NODE_24_length_372190_cov_224.35",
    "NODE_26_length_356957_cov_224.046",
    "NODE_28_length_331107_cov_224.832",
    "NODE_30_length_320616_cov_227.122",
    "NODE_31_length_318595_cov_223.603",
]

# TODO mock ORFs file called by prodigal on test data
# fpath = os.path.join(os.path.dirname(__file__)"/tests/data/test_prots.out"

# fpath = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data/test_prots.out")

# # @patch("test_prodigal.prodigal.get_versions", return_value="2.6")
# def test_contigs_from_headers():
#     with patch("test_prodigal.prodigal.get_versions", return_value="2.0"):
#         out = prodigal.contigs_from_headers(fpath=fpath)
#         assert len(out) == 75390
#         assert type(out) is dict

#     with patch("test_prodigal.prodigal.get_versions", return_value="2.6"):
#         out = prodigal.contigs_from_headers(fpath=fpath)
#         assert len(out) == 75390
#         assert type(out) is dict


# # @patch("test_prodigal.prodigal.get_versions", return_value="2.0")
# def test_orf_records_from_contigs():
#     with patch("test_prodigal.prodigal.get_versions", return_value="2.0"):
#         out = prodigal.orf_records_from_contigs(contigs=contigs, fpath=fpath,)
#         # assert len(out) == 4414

#     with patch("test_prodigal.prodigal.get_versions", return_value="2.6"):
#         out = prodigal.orf_records_from_contigs(contigs=contigs, fpath=fpath,)
#         assert len(out) == 4414


@patch("builtins.open", new_callable=mock_open, read_data="Contents of a fasta file")
# See: https://stackoverflow.com/a/34677735/12671809 it can also be expanded
# as https://queirozf.com/entries/python-unittest-examples-mocking-and-patching#patch-open-file
# @patch("os.path.exists", return_value=True)
@patch("os.path.getsize", return_value=2 * 1024 * 1024)
@patch("autometa.common.external.prodigal.annotate_parallel")
@patch("autometa.common.external.prodigal.annotate_sequential")
def test_output(
    patched_annotate_sequential,
    patched_annotate_parallel,
    patched_file_size,
    patched_open,
):
    nucls_out = MagicMock(
        return_value="path_to_nucs_out", name="nucleic acid_output_file"
    )
    prots_out = MagicMock(
        return_value="path_to_prots_out", name="amino_acid_output_file"
    )
    prodigal.run(
        assembly, nucls_out=nucls_out, prots_out=prots_out, force=True, parallel=False,
    )
    prodigal.run(
        assembly, nucls_out=nucls_out, prots_out=prots_out, force=True, parallel=True,
    )

    assert patched_annotate_parallel.called is True
    assert patched_annotate_sequential.called is True
    assert patched_file_size.called is True
    # assert patched_file_exists.called is True
    assert patched_open.called is True
    assert open("path/to/fasta/file").read() == "Contents of a fasta file"
