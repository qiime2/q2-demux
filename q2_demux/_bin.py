# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import Union
import math

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt
)

single = SingleLanePerSampleSingleEndFastqDirFmt
paired = SingleLanePerSamplePairedEndFastqDirFmt


def bin(
    seqs: Union[single, paired],
    num_bins: int = None,
    bin_delimiters: list[int] = None,
    bin_centers: list[int] = None
) -> dict[str, Union[single, paired]]:

    n_params = bool(num_bins) + bool(bin_delimiters) + bool(bin_centers)
    if n_params == 0:
        raise ValueError(
            "You must pass one of the following parameters: "
            "n_bins, bin_delimiters, or bin_centers."
        )
    elif n_params >= 2:
        raise ValueError(
            "You can't pass more than one of the following parameters: "
            "n_bins, bin_delimiters, or bin_centers."
        )
    if bin_delimiters:
        bins = _parse_bin_delimiters(bin_delimeters=bin_delimiters)
    if bin_centers:
        bins = _parse_bin_centers(bin_centers=bin_centers)


def _parse_bin_centers(bin_centers: list[int]) -> list[tuple]:
    bins = []
    bin_centers.sort()
    for i, center in enumerate(bin_centers):
        bin_delimiter = center + bin_centers[i+1] / 2
        if i == 0:
            bins.append((1, bin_delimiter))
        elif i == len(bin_centers) - 1:
            bins.append((bins[i - 1][1], math.inf))
        else:
            bins.append((bins[i - 1][1], bin_delimiter))
    return bins


def _parse_bin_delimiters(bin_delimeters: list[int]) -> list[tuple]:
    bins = []
    for i, delim in enumerate(bin_delimeters):
        if i == 0:
            bins.append((1, delim))
        elif i == len(bin_delimeters) - 1:
            bins.append((delim, math.inf))
        else:
            bins.append((bins[i - 1][1], delim))


def _parse_num_bins():
    pass
