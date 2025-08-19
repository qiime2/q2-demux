# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import gzip
import os

import pandas as pd

from qiime2 import Metadata
from qiime2.util import duplicate
from q2_types.per_sample_sequences import \
    CasavaOneEightSingleLanePerSampleDirFmt

from ._summarize import _PlotQualView


def filter_samples(demux: _PlotQualView, metadata: Metadata = None,
                   where: str = None, exclude_ids: bool = False,
                   remove_empty: bool = False) \
                   -> CasavaOneEightSingleLanePerSampleDirFmt:

    if not any([metadata, remove_empty]):
        raise ValueError(
            "At least one of the following parameters must be provided: "
            "metadata, remove-empty."
        )

    results = CasavaOneEightSingleLanePerSampleDirFmt()
    paired = demux.paired
    samples = demux.directory_format
    manifest = samples.manifest.view(pd.DataFrame)
    ids_to_keep = set(manifest.index)

    if metadata is not None:
        ids_to_keep = metadata.get_ids(where=where)
        if not ids_to_keep:
            raise ValueError('No filtering requested.')

    if exclude_ids:
        ids_to_keep = set(manifest.index) - set(ids_to_keep)

    if remove_empty:
        ids_empty = _get_empty_sample_ids(manifest)
        ids_to_keep = set(ids_to_keep) - set(ids_empty)

    try:
        for id in ids_to_keep:
            forward = manifest.loc[id].forward
            duplicate(forward, os.path.join(str(results),
                      os.path.split(forward)[1]))
            if paired:
                reverse = manifest.loc[id].reverse
                duplicate(reverse, os.path.join(str(results),
                          os.path.split(reverse)[1]))
    except KeyError:
        raise ValueError(f'{id!r} is not a sample present in the '
                         'demultiplexed data.')

    return results


def _get_empty_sample_ids(manifest):
    """
    Identify and return sample names from a manifest DataFrame where
    at least one FASTQ file is empty.

    Parameters:
        manifest (pandas.DataFrame): A DataFrame where the index contains
        sample names and the values are file paths pointing to FASTQ files.

    Returns:
        list: A list of sample names for which at least one associated
        file is empty (no content on the first line).
    """
    empty_samples = []
    for sample, row in manifest.iterrows():
        for path in row:
            if path is not None:
                with gzip.open(path, "rt") as f:
                    if not f.readline():
                        empty_samples.append(sample)
                        break
    return empty_samples
