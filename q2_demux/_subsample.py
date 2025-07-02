# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
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

from ._util import read_fastq_seqs


def subsample_single(
                     sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                     fraction: float,
                     drop_empty: bool = False
                     ) -> CasavaOneEightSingleLanePerSampleDirFmt:
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest.view(pd.DataFrame)
    for _, fwd_path in manifest.itertuples():
        fwd_name = os.path.basename(fwd_path)
        fwd_path_in = str(sequences.path / fwd_name)
        fwd_path_out = str(result.path / fwd_name)

        with gzip.open(str(fwd_path_out), mode='w') as fwd:
            for fwd_rec in read_fastq_seqs(fwd_path_in):
                if random.random() <= fraction:
                    fwd.write(('\n'.join(fwd_rec) + '\n').encode('utf-8'))

    if drop_empty:
        remove_empty_files(
            sequences,
            result
        )

    return result


def subsample_paired(
                     sequences: SingleLanePerSamplePairedEndFastqDirFmt,
                     fraction: float,
                     drop_empty: bool = False
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
                file_pair = zip(read_fastq_seqs(fwd_path_in),
                                read_fastq_seqs(rev_path_in))
                for fwd_rec, rev_rec in file_pair:
                    if random.random() <= fraction:
                        fwd.write(
                            ('\n'.join(fwd_rec) + '\n').encode('utf-8'))
                        rev.write(
                            ('\n'.join(rev_rec) + '\n').encode('utf-8'))

    if drop_empty:
        remove_empty_files(
            sequences,
            result
        )

    return result


def remove_empty_files(sequences: SingleLanePerSamplePairedEndFastqDirFmt,
                       result: CasavaOneEightSingleLanePerSampleDirFmt
                       ):
    """
    This function removes files from the `result` directory if there are no
    reads after random subsampling. Afterward the files are also removed
    from the MANIFEST file.
    """
    empty_files = []

    file_list = os.listdir(str(result))
    sf_path = result.path
    mf_path_in = str(sequences.path / 'MANIFEST')
    mf_path_out = str(result.path / 'MANIFEST')

    for file in file_list:
        file_path = sf_path / file
        gz_file = gzip.GzipFile(str(file_path), 'rb')
        if gz_file.peek(1) == b'':
            os.remove(sf_path / file)
            empty_files.append(file)

    with open(mf_path_in, mode='r') as mf:
        lines = mf.readlines()

    new_lines = []
    for line in lines:
        for empty_file in empty_files:
            if empty_file in line:
                break
        else:
            new_lines.append(line)

    with open(mf_path_out, mode='w') as mf:
        mf.writelines(new_lines)

    if len(os.listdir(result.path)) == 1:
        raise ValueError('All sample were empty after subsampling, try again'
                         ' with a larger fraction')
