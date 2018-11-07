# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import gzip
import random

import pandas as pd

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt)

from q2_demux._demux import _read_fastq_seqs


def subsample_single(sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                     fraction: float
                     ) -> CasavaOneEightSingleLanePerSampleDirFmt:
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest.view(pd.DataFrame)

    for _, fwd_path in manifest.itertuples():
        fwd_name = os.path.basename(fwd_path)
        fwd_path_in = str(sequences.path / fwd_name)
        fwd_path_out = str(result.path / fwd_name)
        with gzip.open(str(fwd_path_out), mode='w') as fwd:
            for fwd_rec in _read_fastq_seqs(fwd_path_in):
                if random.random() <= fraction:
                    fwd.write(('\n'.join(fwd_rec) + '\n').encode('utf-8'))

    return result


def subsample_paired(sequences: SingleLanePerSamplePairedEndFastqDirFmt,
                     fraction: float
                     ) -> CasavaOneEightSingleLanePerSampleDirFmt:
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest.view(pd.DataFrame)

    for _, fwd_path, rev_path in manifest.itertuples():
        fwd_name = os.path.basename(fwd_path)
        rev_name = os.path.basename(rev_path)
        fwd_path_in = str(sequences.path / fwd_name)
        rev_path_in = str(sequences.path / rev_name)
        fwd_path_out = str(result.path / fwd_name)
        rev_path_out = str(result.path / rev_name)
        with gzip.open(str(fwd_path_out), mode='w') as fwd:
            with gzip.open(str(rev_path_out), mode='w') as rev:
                file_pair = zip(_read_fastq_seqs(fwd_path_in),
                                _read_fastq_seqs(rev_path_in))
                for fwd_rec, rev_rec in file_pair:
                    if random.random() <= fraction:
                        fwd.write(('\n'.join(fwd_rec) + '\n').encode('utf-8'))
                        rev.write(('\n'.join(rev_rec) + '\n').encode('utf-8'))

    return result
