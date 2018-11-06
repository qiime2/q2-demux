# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import gzip
import yaml
import random

import pandas as pd

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqManifestFormat, YamlFormat)

# def _subsample_paired(fastq_map):
#     qual_sample = collections.defaultdict(list)
#     min_seq_len = {'forward': float('inf'), 'reverse': float('inf')}
#     for fwd, rev, index in fastq_map:
#         file_pair = zip(_read_fastq_seqs(fwd), _read_fastq_seqs(rev))
#         for i, (fseq, rseq) in enumerate(file_pair):
#             if i == index[0]:
#                 min_seq_len['forward'] = min(min_seq_len['forward'],
#                                              len(fseq[1]))
#                 min_seq_len['reverse'] = min(min_seq_len['reverse'],
#                                              len(rseq[1]))
#                 qual_sample['forward'].append(_decode_qual_to_phred33(fseq[3]))
#                 qual_sample['reverse'].append(_decode_qual_to_phred33(rseq[3]))
#                 index.pop(0)
#                 if len(index) == 0:
#                     break
#     return qual_sample, min_seq_len
#
#
# def _subsample_single(fastq_map):
#     qual_sample = collections.defaultdict(list)
#     min_seq_len = {'forward': float('inf'), 'reverse': None}
#     for file, index in fastq_map:
#         for i, seq in enumerate(_read_fastq_seqs(file)):
#             if i == index[0]:
#                 min_seq_len['forward'] = min(min_seq_len['forward'],
#                                              len(seq[1]))
#                 qual_sample['forward'].append(_decode_qual_to_phred33(seq[3]))
#                 index.pop(0)
#                 if len(index) == 0:
#                     break
#     return qual_sample, min_seq_len

def subsample_single(sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                     percentage: float
                    ) -> SingleLanePerSampleSingleEndFastqDirFmt:
    # result = SingleLanePerSamplePairedEndFastqDirFmt()
    #
    # manifest = FastqManifestFormat()
    # manifest_fh = manifest.open()
    # manifest_fh.write('sample-id,filename,direction\n')

    print(dir(sequences.sequences))


    return sequences

def subsample_paired(sequences: SingleLanePerSamplePairedEndFastqDirFmt,
                     percentage: float
                    ) -> SingleLanePerSamplePairedEndFastqDirFmt:
    result = SingleLanePerSamplePairedEndFastqDirFmt()

    manifest_path = os.path.join(str(sequences),sequences.manifest.pathspec)

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
                fwd.write(gzip.open(fwd_path_in).read())
                rev.write(gzip.open(rev_path_in).read())

    _write_metadata_yaml(result)

    manifest = FastqManifestFormat()
    with manifest.open() as manifest_fh:
        manifest_fh.write(open(manifest_path).read())
    result.manifest.write_data(manifest, FastqManifestFormat)

    return result

def _write_metadata_yaml(dir_fmt):
    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    dir_fmt.metadata.write_data(metadata, YamlFormat)
