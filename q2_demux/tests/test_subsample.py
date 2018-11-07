# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import itertools
import gzip
import unittest

import pandas as pd
import numpy as np

import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.per_sample_sequences import (
    FastqGzFormat, CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_demux._demux import (BarcodePairedSequenceFastqIterator)
from q2_demux import (emp_paired,
                      subsample_single, subsample_paired)


class SubsampleTests(TestPluginBase):
    # this functionality is derived from test_demux.EmpTestingUtils
    package = 'q2_demux.tests'

    def _get_total_sequence_count(self, seq_ids):
        return len(list(itertools.chain(*seq_ids)))

    def _validate_fastq_subsampled(self, obs, exp, forward : bool):
        subsampled_sequence_ids = []
        observed_samples = 0

        obs = obs.sequences.iter_views(FastqGzFormat)
        exp = exp.sequences.iter_views(FastqGzFormat)

        # Iterate over each sample, side-by-side
        for (_, exp_fp), (_, obs_fp) in zip(exp, obs):
            # Process fwd only if `forward==True`
            if forward and 'R2' in str(exp_fp):
                continue
            # Process rev only if `forward==False`
            if not forward and 'R1' in str(exp_fp):
                continue

            observed_samples += 1
            exp_fh = gzip.open(str(exp_fp), 'rt')
            obs_fh = gzip.open(str(obs_fp), 'rt')

            # Assemble expected sequences, per-sample
            exp_seqs = [r for r in itertools.zip_longest(*[exp_fh] * 4)]

            # Assemble observed sequences, per-sample
            obs_seqs = [r for r in itertools.zip_longest(*[obs_fh] * 4)]

            # the number of output sequences is less than or equal to the
            # number of input sequences
            self.assertTrue(len(obs_seqs) <= len(exp_seqs))

            # is the observed set a subset of expected?
            self.assertTrue(set(obs_seqs).issubset(set(exp_seqs)))

            subsampled_sequence_ids.append([e[0] for e in obs_seqs])

        # return the output sequence IDs, so that they can be used in
        # other tests
        return subsampled_sequence_ids, observed_samples


class SubsampleSingleTests(SubsampleTests):

    def setUp(self):
        super().setUp()

        demuxed = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path('subsample_single_end'), mode='r')
        self.demux_data = transform(
            demuxed, to_type=SingleLanePerSampleSingleEndFastqDirFmt)

    def test_subsample(self):
        actual = subsample_single(self.demux_data, fraction=0.5)

        fwd_subsampled_sequence_ids, obs_sample_count = \
            self._validate_fastq_subsampled(actual, self.demux_data,
                                            forward=True)

        self.assertEqual(obs_sample_count, 5)

        # some sequences have been removed - this could occasionally fail,
        # but the frequency of that should be ~ 2 * 0.5 ** 11
        seq_count = self._get_total_sequence_count(fwd_subsampled_sequence_ids)
        self.assertTrue(0 < seq_count < 11)

    def test_correct_output_files_on_small_subsample(self):
        # some or all of the output files are likely to be empty, but they
        # should still be present and in the manifest
        actual = subsample_single(self.demux_data, fraction=0.00001)

        _, obs_sample_count = self._validate_fastq_subsampled(actual,
                                                              self.demux_data,
                                                              forward=True)

        self.assertEqual(obs_sample_count, 5)


class SubsamplePairedTests(SubsampleTests):

    def setUp(self):
        super().setUp()

        demuxed = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path('subsample_paired_end'), mode='r')
        self.demux_data = transform(
            demuxed, to_type=SingleLanePerSamplePairedEndFastqDirFmt)

    def test_subsample(self):
        actual = subsample_paired(self.demux_data, fraction=0.5)

        fwd_subsampled_sequence_ids, fwd_obs_sample_count = \
            self._validate_fastq_subsampled(actual, self.demux_data,
                                            forward=True)
        rev_subsampled_sequence_ids, rev_obs_sample_count = \
            self._validate_fastq_subsampled(actual, self.demux_data,
                                            forward=False)

        self.assertEqual(fwd_obs_sample_count, 5)
        self.assertEqual(rev_obs_sample_count, 5)

        # some sequences have been removed - this could occasionally fail,
        # but the frequency of that should be ~ 2 * 0.5 ** 11
        f_seq_count = self._get_total_sequence_count(fwd_subsampled_sequence_ids)
        self.assertTrue(0 < f_seq_count < 11)
        r_seq_count = self._get_total_sequence_count(rev_subsampled_sequence_ids)
        self.assertTrue(0 < r_seq_count < 11)
        self.assertEqual(f_seq_count, r_seq_count)

        print(fwd_subsampled_sequence_ids)
        print(rev_subsampled_sequence_ids)
        self.assertEqual(fwd_subsampled_sequence_ids,
                         rev_subsampled_sequence_ids)

    def test_correct_output_files_on_small_subsample(self):
        # some or all of the output files are likely to be empty, but they
        # should still be present and in the manifest
        actual = subsample_paired(self.demux_data, fraction=0.00001)

        _, fwd_obs_sample_count = \
            self._validate_fastq_subsampled(actual, self.demux_data,
                                            forward=True)
        _, rev_obs_sample_count = \
            self._validate_fastq_subsampled(actual, self.demux_data,
                                            forward=False)

        self.assertEqual(fwd_obs_sample_count, 5)
        self.assertEqual(rev_obs_sample_count, 5)


if __name__ == '__main__':
    unittest.main()
