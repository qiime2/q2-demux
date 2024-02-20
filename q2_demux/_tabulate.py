# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import pandas as pd

import qiime2
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt)

from q2_types.feature_data._util import _read_fastq_seqs


def tabulate_read_counts(sequences: SingleLanePerSampleSingleEndFastqDirFmt
                         ) -> qiime2.Metadata:
    result = {}

    for e in sequences:
        manifest = e.manifest.view(pd.DataFrame)
        for record in manifest.itertuples():
            sample_id = record[0]
            fwd_path = record[1]
            read_count = 0
            if sample_id in result:
                raise KeyError("At least one duplicated sample id was "
                               f"detected ({sample_id}). "
                               "Sample ids must be unique across inputs.")
            fwd_name = os.path.basename(fwd_path)
            fwd_path = str(e.path / fwd_name)
            for fwd_rec in _read_fastq_seqs(fwd_path):
                read_count += 1
            result[sample_id] = read_count

    result = pd.Series(result)
    result.name = 'Demultiplexed sequence count'
    result = result.to_frame()
    result.index.name = 'sample-id'

    return qiime2.Metadata(result)
