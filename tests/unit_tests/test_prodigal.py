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


def test_output():
    nucs_out, prots_out = prodigal.run(
        "/home/the_bio_informatician/Autometa/tests/data/metagenome.fna",
        "nucs.out",
        "prots.out",
    )
    assert nucs_out == "nucs.out"
    assert prots_out == "prots.out"
    assert type(nucs_out) is str
    assert type(prots_out) is str
