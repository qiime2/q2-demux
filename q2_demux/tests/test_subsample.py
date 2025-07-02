# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import itertools
import gzip
import unittest
import os

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.per_sample_sequences import (
    FastqGzFormat, CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_demux import subsample_single, subsample_paired


class SubsampleTests(TestPluginBase):
    # this functionality is derived from test_demux.EmpTestingUtils
    package = 'q2_demux.tests'

    def _get_total_sequence_count(self, seq_ids):
        return len(list(itertools.chain(*seq_ids)))

    def _validate_fastq_subsampled(self, obs, exp, forward):
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

    def test_subsample_single_drop_empty_reads_all(self):
        """
        This function tests that if the subsampling fraction is zero an error
        is raised alerting the user that all samples are empty.
        """
        path = self.get_data_path('subsample_data_test_single')
        view = SingleLanePerSampleSingleEndFastqDirFmt(path, mode='r')

        with self.assertRaisesRegex(
            ValueError,
            'All sample were empty after subsampling, try again with a larger '
            'fraction.'
        ):
            subsample_single(view, fraction=0.0, drop_empty=True)

    def test_subsample_single_drop_empty_reads_some(self):
        """
        This function tests that if some samples are empty and some are not
        after subsampling, only the empty samples are dropped for single end
        reads.
        """
        path = self.get_data_path('subsample_data_test_single')
        view = SingleLanePerSampleSingleEndFastqDirFmt(path, mode='r')

        dropped_some = subsample_single(
            view, fraction=0.0008, drop_empty=True
        )

        path = dropped_some.path
        file_path_removed = 'sample-short_S2_L001_R1_001.fastq.gz'
        path_removed = os.path.join(path, file_path_removed)
        file_path_kept = 'sample-long_S1_L001_R1_001.fastq.gz'
        path_kept = os.path.join(path, file_path_kept)

        try:
            self.assertFalse(os.path.exists(path_removed))
            self.assertTrue(os.path.exists(path_kept))

        except AssertionError:
            raise AssertionError(
                "This test fails approximately 1 in 1000 times. Run the test "
                "again."
            )


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
        f_seq_count = self._get_total_sequence_count(
            fwd_subsampled_sequence_ids)
        r_seq_count = self._get_total_sequence_count(
            rev_subsampled_sequence_ids)
        self.assertTrue(0 < f_seq_count < 11)
        self.assertTrue(0 < r_seq_count < 11)

        self.assertEqual(f_seq_count, r_seq_count)

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

    def test_subsample_paired_drop_empty_reads_all(self):
        """
        This function tests that if the subsampling fraction is zero an error
        will be raised alerting the user that all samples are empty.
        """
        path = self.get_data_path('subsample_data_test_paired')
        view = SingleLanePerSamplePairedEndFastqDirFmt(path, mode='r')

        with self.assertRaisesRegex(
            ValueError,
            'All sample were empty after subsampling, try again with a larger '
            'fraction.'
        ):
            subsample_paired(view, fraction=0.0, drop_empty=True)

    def test_subsample_paired_drop_empty_reads_some(self):
        """
        This function tests that if some samples are empty and some are not
        after subsampling, only the empty samples are dropped for paired end
        reads.
        """
        path = self.get_data_path('subsample_data_test_paired')
        view = SingleLanePerSamplePairedEndFastqDirFmt(path, mode='r')

        dropped_some = subsample_paired(
            view, fraction=0.0008, drop_empty=True
        )

        path = dropped_some.path
        file_path = 'sample2_2_L001_R1_001.fastq.gz'
        path_remove = os.path.join(path, file_path)
        file_path_rev = 'sample2_2_L001_R2_001.fastq.gz'
        path_remove_rev = os.path.join(path, file_path_rev)

        file_path_keep = 'sample1_1_L001_R1_001.fastq.gz'
        path_keep = os.path.join(path, file_path_keep)
        file_path_keep_rev = 'sample1_1_L001_R2_001.fastq.gz'
        path_keep_rev = os.path.join(path, file_path_keep_rev)

        try:
            self.assertFalse(os.path.exists(path_remove))
            self.assertFalse(os.path.exists(path_remove_rev))
            self.assertTrue(os.path.exists(path_keep))
            self.assertTrue(os.path.exists(path_keep_rev))

        except AssertionError:
            raise AssertionError("This test fails approximately 1 in 1000 "
                                 "times. Run the test again")


if __name__ == '__main__':
    unittest.main()
