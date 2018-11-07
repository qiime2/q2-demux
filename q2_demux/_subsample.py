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
    CasavaOneEightSingleLanePerSampleDirFmt,
    FastqManifestFormat)

from q2_demux._demux import _read_fastq_seqs, _write_metadata_yaml


def subsample_single(sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                     fraction: float
                     ) -> CasavaOneEightSingleLanePerSampleDirFmt:
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest.view(pd.DataFrame)

    for _, fwd_name in manifest.itertuples():
        # can the following be modified?
        fwd_path_in = os.path.join(str(sequences), fwd_name)
        fwd_path_out = str(result.path / os.path.basename(fwd_name))
        with gzip.open(str(fwd_path_out), mode='w') as fwd:
            for fwd_rec in _read_fastq_seqs(fwd_path_in):
                if random.random() <= fraction:
                    fwd.write(('\n'.join(fwd_rec) + '\n').encode('utf-8'))

    return result


def subsample_paired(sequences: SingleLanePerSamplePairedEndFastqDirFmt,
                     fraction: float
                     ) -> SingleLanePerSamplePairedEndFastqDirFmt:
    result = SingleLanePerSamplePairedEndFastqDirFmt()

    manifest_path = os.path.join(str(sequences), sequences.manifest.pathspec)

    manifest = pd.read_csv(manifest_path, header=0, comment='#')

    fwd = manifest[manifest.direction == 'forward'].filename.tolist()
    rev = manifest[manifest.direction == 'reverse'].filename.tolist()
    sample_map = [(file, rev[fwd.index(file)]) for file in fwd]

    for fwd_name, rev_name in sample_map:
        fwd_path_in = os.path.join(str(sequences), fwd_name)
        rev_path_in = os.path.join(str(sequences), rev_name)
        sample_id, barcode_id, _ = fwd_name.split('_', maxsplit=2)
        fwd_path_out = result.sequences.path_maker(sample_id=sample_id,
                                                   barcode_id=barcode_id,
                                                   lane_number=1,
                                                   read_number=1)
        rev_path_out = result.sequences.path_maker(sample_id=sample_id,
                                                   barcode_id=barcode_id,
                                                   lane_number=1,
                                                   read_number=2)
        with gzip.open(str(fwd_path_out), mode='w') as fwd:
            with gzip.open(str(rev_path_out), mode='w') as rev:
                file_pair = zip(_read_fastq_seqs(fwd_path_in),
                                _read_fastq_seqs(rev_path_in))
                for fwd_rec, rev_rec in file_pair:
                    if random.random() <= fraction:
                        fwd.write(('\n'.join(fwd_rec) + '\n').encode('utf-8'))
                        rev.write(('\n'.join(rev_rec) + '\n').encode('utf-8'))

    _write_metadata_yaml(result)

    manifest = FastqManifestFormat()
    with manifest.open() as manifest_fh:
        manifest_fh.write(open(manifest_path).read())
    result.manifest.write_data(manifest, FastqManifestFormat)

    return result
