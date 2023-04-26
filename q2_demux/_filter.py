# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import pandas as pd

from qiime2 import Metadata
from qiime2.util import duplicate
from q2_types.per_sample_sequences import \
    CasavaOneEightSingleLanePerSampleDirFmt

from ._summarize import _PlotQualView


def filter_samples(demux: _PlotQualView, metadata: Metadata,
                   where: str = None, exclude_ids: bool = False) \
                   -> CasavaOneEightSingleLanePerSampleDirFmt:
    results = CasavaOneEightSingleLanePerSampleDirFmt()

    paired = demux.paired
    samples = demux.directory_format

    ids_to_keep = metadata.get_ids(where=where)
    if not ids_to_keep:
        raise ValueError('No filtering requested.')
    manifest = samples.manifest.view(pd.DataFrame)
    # The ids we read out of the metadata are strings, so we need to make sure
    # the manifest is indexed on strings. If we didn't do this, you could not
    # use numbers as ids because they would be read from the metadata as
    # strings, but the manifest's index would be numeric. This would lead to
    # errors when we do manifest.loc[id] because we would be using a string
    # index into a numeric dataframe.
    manifest.index = manifest.index.astype(str)

    if exclude_ids:
        ids_to_keep = set(manifest.index) - set(ids_to_keep)

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
