#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pytest

from autometa.common.external import prodigal

assembly = "tests/data/metagenome.fna"


# @pytest.fixture(scope="session")
# def test_prodigal(tmpdir_factory):
#     tmpdir = tmpdir_factory.mktemp("output")
#     nucls_out = os.path.join(tmpdir, "nucls.out")
#     prots_out = os.path.join(tmpdir, "prots.out")
#     with pytest.raises(FileNotFoundError):
#         prodigal.run(
#             assembly=assembly, nucls_out=nucls_out, prots_out=prots_out,
#         )

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
fpath = "/home/the_bio_informatician/Autometa/tests/data/test_prots.out"
# TODO mock prodigal version
def test_contigs_from_headers():
    out = prodigal.contigs_from_headers(fpath=fpath)
    assert len(out) == 75390
    assert type(out) is dict


def test_orf_records_from_contigs():
    out = prodigal.orf_records_from_contigs(contigs=contigs, fpath=fpath,)
    assert len(out) == 4414


# def test_output():
#     nucs_out, prots_out = prodigal.run(
#         "/home/the_bio_informatician/Autometa/tests/data/metagenome.fna",
#         "nucs.out",
#         "prots.out",
#     )
#     assert nucs_out == "nucs.out"
#     assert prots_out == "prots.out"
#     assert type(nucs_out) is str
#     assert type(prots_out) is str
