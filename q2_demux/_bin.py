# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import Union
import math
import pandas as pd
import os


from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt
)

from ._util import read_fastq_seqs

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
    elif bin_centers:
        bins = _parse_bin_centers(bin_centers=bin_centers)
    elif num_bins:
        bins = _parse_num_bins(seqs=seqs, num_bins=num_bins)




def _parse_bin_centers(bin_centers: list[int]) -> list[tuple]:
    bins = []
    bin_centers.sort()
    for i, center in enumerate(bin_centers):
        if i == len(bin_centers) - 1:
            bins.append((bins[i - 1][1], math.inf))
            continue
        bin_delimiter = center + bin_centers[i+1] / 2
        if i == 0:
            bins.append((1, bin_delimiter))
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
    return bins


def _parse_num_bins(seqs: Union[single, paired], num_bins: int) -> list[tuple]:
    bins = []
    base_fp = str(seqs)
    manifest_df = seqs.manifest.view(pd.DataFrame)
    seq_max = 0
    seq_min = math.inf

    for _, fp in manifest_df.iterrows():
        ffp = fp["forward"]
        rfp = None
        for _, sequence, _, _ in read_fastq_seqs(os.path.join(base_fp, ffp)):
            if len(sequence) < seq_min:
                seq_min = len(sequence)
            if len(sequence) > seq_max:
                seq_max = len(sequence)

        if "reverse" in fp:
            rfp = fp["reverse"]
            for _, sequence, _, _ in read_fastq_seqs(os.path.join(base_fp,
                                                                  rfp)):
                if len(sequence) < seq_min:
                    seq_min = len(sequence)
                if len(sequence) > seq_max:
                    seq_max = len(sequence)
    seq_range = seq_max - seq_min
    bin_range = seq_range//num_bins
    for i in num_bins:
        if i == 0:
            bins.append((seq_min, seq_min + bin_range))
        elif i == len(num_bins) - 1:
            bins.append((bins[i - 1][1], math.inf))
        else:
            bins.append((bins[i - 1][1], bins[i - 1][1] + bin_range))
    return bins
