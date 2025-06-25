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


def subsample_single(sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                     fraction: float,
                     drop_empty: bool = False
                     ) -> CasavaOneEightSingleLanePerSampleDirFmt:
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest.view(pd.DataFrame)
    if drop_empty:
        result = remove_files(sequences, fraction)
    else:
        for _, fwd_path in manifest.itertuples():
            fwd_name = os.path.basename(fwd_path)
            fwd_path_in = str(sequences.path / fwd_name)
            fwd_path_out = str(result.path / fwd_name)

            with gzip.open(str(fwd_path_out), mode='w') as fwd:
                for fwd_rec in read_fastq_seqs(fwd_path_in):
                    if random.random() <= fraction:
                        fwd.write(('\n'.join(fwd_rec) + '\n').encode('utf-8'))

    return result


def subsample_paired(sequences: SingleLanePerSamplePairedEndFastqDirFmt,
                     fraction: float,
                     drop_empty: bool = False
                     ) -> CasavaOneEightSingleLanePerSampleDirFmt:
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest.view(pd.DataFrame)

    if drop_empty:
        result = remove_files(sequences, fraction, is_paired=True)
    else:
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

    return result


def remove_files(sequences,
                 fraction: float,
                 is_paired: bool = False,
                 ) -> CasavaOneEightSingleLanePerSampleDirFmt:
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest.view(pd.DataFrame)
    removed_files = []
    mf_path_in = str(sequences.path / 'MANIFEST')
    mf_path_out = str(result.path / 'MANIFEST')

    if is_paired:
        for _, fwd_path, rev_path in manifest.itertuples():
            fwd_name = os.path.basename(fwd_path)
            rev_name = os.path.basename(rev_path)
            fwd_path_in = str(sequences.path / fwd_name)
            rev_path_in = str(sequences.path / rev_name)
            fwd_path_out = str(result.path / fwd_name)
            rev_path_out = str(result.path / rev_name)
            reads = 0

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
                            reads += 1
            if reads == 0:
                os.remove(fwd_path_out)
                os.remove(rev_path_out)
                removed_files.append(fwd_name)
                removed_files.append(rev_name)
    else:
        for _, fwd_path in manifest.itertuples():
            fwd_name = os.path.basename(fwd_path)
            fwd_path_in = str(sequences.path / fwd_name)
            fwd_path_out = str(result.path / fwd_name)

            reads = 0
            with gzip.open(str(fwd_path_out), mode='w') as fwd:
                for fwd_rec in read_fastq_seqs(fwd_path_in):
                    if random.random() <= fraction:
                        fwd.write(('\n'.join(fwd_rec) + '\n').encode('utf-8'))
                        reads += 1
            if reads == 0:
                os.remove(fwd_path_out)
                removed_files.append(fwd_name)

    with open(mf_path_in, mode='r') as mf:
        lines = mf.readlines()

    new_lines = []
    for line in lines:
        keep_file = True
        for removed_file in removed_files:
            if removed_file in line:
                keep_file = False
                break
        if keep_file:
            new_lines.append(line)

    with open(mf_path_out, mode='w') as mf:
        mf.writelines(new_lines)

    return result
