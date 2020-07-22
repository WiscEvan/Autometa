#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
from autometa.common.external import samtools

# TODO: Add mock for sam file
sam_fpath = "<path/to/mock/sam/file>"


def test_sort():

    with pytest.raises(FileNotFoundError):
        samtools.sort(sam="sam", bam="bam")
    with pytest.raises(TypeError):
        samtools.sort(sam="sam_fpath", bam="bam", nproc=2.9)
    with pytest.raises(TypeError):
        samtools.sort(
            sam="sam_fpath", bam="bam", nproc=-2,
        )
