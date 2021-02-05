# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import os

import pandas as pd

from qiime2 import Metadata
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)

from q2_demux._filter import filter_samples
from q2_demux._summarize import _PlotQualView


class FilterSamplesTests(TestPluginBase):
    package = 'q2_demux.tests'

    def setUp(self):
        super().setUp()

        data_single = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('filter_samples_single_end/dir_fmt'), mode='r')
        self.sample_single = _PlotQualView(data_single, False)
        self.manifest_single = data_single.manifest.view(pd.DataFrame)

        self.md_single_all = Metadata.load(
            self.get_data_path('filter_samples_single_end/filter_all.tsv'))
        self.md_single_subset = Metadata.load(
            self.get_data_path('filter_samples_single_end/filter_subset.tsv'))
        self.md_single_none = Metadata.load(
            self.get_data_path('filter_samples_single_end/filter_none.tsv'))

        data_paired = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path('filter_samples_paired_end/dir_fmt'), mode='r')
        self.sample_paired = _PlotQualView(data_paired, True)
        self.manifest_paired = data_paired.manifest.view(pd.DataFrame)

        self.md_paired_all = Metadata.load(
            self.get_data_path('filter_samples_single_end/filter_all.tsv'))
        self.md_paired_subset = Metadata.load(
            self.get_data_path('filter_samples_single_end/filter_subset.tsv'))
        self.md_paired_none = Metadata.load(
            self.get_data_path('filter_samples_single_end/filter_none.tsv'))

    def _assert_single_contains(self, dir_fmt, exp_ids):
        obs_ids = [file.split('_', 1)[0] for file in
                   os.listdir(str(dir_fmt))]
        self.assertEqual(set(obs_ids), set(exp_ids))

    def _assert_paired_contains(self, dir_fmt, exp_ids):
        obs_ids = [(file.split('_')[0::3]) for file in
                   os.listdir(str(dir_fmt))]
        obs_ids = [(id[0] + id[1]) for id in obs_ids]
        self.assertEqual(set(obs_ids), set(exp_ids))

    def test_filter_single_all(self):
        exps = [(None, ['sample1', 'sample2']),
                ("Study='A'", ['sample1']),
                ("Study='B'", ['sample2']),
                ("Study='A' OR Study='B'", ['sample1', 'sample2'])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_single, self.md_single_all,
                                     where, False)
            self._assert_single_contains(dir_fmt, exp)

    def test_filter_single_all_exclude(self):
        exps = [(None, []),
                ("Study='A'", ['sample2']),
                ("Study='B'", ['sample1']),
                ("Study='A' OR Study='B'", [])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_single, self.md_single_all,
                                     where, True)
            self._assert_single_contains(dir_fmt, exp)

    def test_filter_single_subset(self):
        exps = [(None, ['sample1']),
                ("Study='A'", ['sample1']),
                ("Study='A' OR Study='B'", ['sample1'])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_single,
                                     self.md_single_subset, where, False)
            self._assert_single_contains(dir_fmt, exp)

    def test_filter_single_subset_no_filter(self):
        wheres = ["Study='B'", "Study='A' AND Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, 'No'):
                filter_samples(self.sample_single, self.md_single_subset,
                               where, False)

    def test_filter_single_subset_exclude(self):
        exps = [(None, ['sample2']),
                ("Study='A'", ['sample2']),
                ("Study='A' OR Study='B'", ['sample2'])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_single,
                                     self.md_single_subset, where, True)
            self._assert_single_contains(dir_fmt, exp)

    def test_filter_sinlge_subset_exclude_no_filter(self):
        wheres = ["Study='B'", "Study='A' AND Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, 'No'):
                filter_samples(self.sample_single, self.md_single_subset,
                               where, True)

    def test_filter_single_none_exclude(self):
        exps = [(None, ['sample1', 'sample2']),
                ("Study='A'", ['sample1', 'sample2']),
                ("Study='A' OR Study='B'", ['sample1', 'sample2'])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_single, self.md_single_none,
                                     where, True)
            self._assert_single_contains(dir_fmt, exp)

    def test_filter_single_id_not_present(self):
        wheres = [None, "Study='A'", "Study='A' OR Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, '1'):
                filter_samples(self.sample_single, self.md_single_none, where,
                               False)

    def test_filter_single_none_no_filter(self):
        wheres = ["Study='B'", "Study='A' AND Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, 'No'):
                filter_samples(self.sample_single, self.md_single_none, where,
                               False)

    def test_filter_single_none_exclude_no_filter(self):
        wheres = ["Study='B'", "Study='A' AND Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, 'No'):
                filter_samples(self.sample_single, self.md_single_none, where,
                               True)

    def test_filter_paired_all(self):
        exps = [(None, ['sample1R1', 'sample1R2', 'sample2R1',
                'sample2R2']),
                ("Study='A'", ['sample1R1', 'sample1R2']),
                ("Study='B'", ['sample2R1', 'sample2R2']),
                ("Study='A' OR Study='B'", ['sample1R1', 'sample1R2',
                 'sample2R1', 'sample2R2'])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_paired, self.md_paired_all,
                                     where, False)
            self._assert_paired_contains(dir_fmt, exp)

    def test_filter_paired_all_exclude(self):
        exps = [(None, []),
                ("Study='A'", ['sample2R1', 'sample2R2']),
                ("Study='B'", ['sample1R1', 'sample1R2']),
                ("Study='A' OR Study='B'", [])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_paired, self.md_paired_all,
                                     where, True)
            self._assert_paired_contains(dir_fmt, exp)

    def test_filter_paired_subset(self):
        exps = [(None, ['sample1R1', 'sample1R2']),
                ("Study='A'", ['sample1R1', 'sample1R2']),
                ("Study='A' OR Study='B'", ['sample1R1', 'sample1R2'])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_paired,
                                     self.md_paired_subset, where, False)
            self._assert_paired_contains(dir_fmt, exp)

    def test_filter_paired_subset_no_filtering(self):
        wheres = ["Study='B'", "Study='A' AND Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, 'No'):
                filter_samples(self.sample_paired, self.md_paired_subset,
                               where, False)

    def test_filter_paired_subset_exclude(self):
        exps = [(None, ['sample2R1', 'sample2R2']),
                ("Study='A'", ['sample2R1', 'sample2R2']),
                ("Study='A' OR Study='B'", ['sample2R1', 'sample2R2'])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_paired,
                                     self.md_paired_subset, where, True)
            self._assert_paired_contains(dir_fmt, exp)

    def test_filter_paired_subset_exclude_no_filtering(self):
        wheres = ["Study='B'", "Study='A' AND Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, 'No'):
                filter_samples(self.sample_paired, self.md_paired_subset,
                               where, True)

    def test_filter_paired_none_exclude(self):
        exps = [(None, ['sample1R1', 'sample1R2', 'sample2R1', 'sample2R2']),
                ("Study='A'", ['sample1R1', 'sample1R2', 'sample2R1',
                               'sample2R2']),
                ("Study='A' OR Study='B'", ['sample1R1', 'sample1R2',
                                            'sample2R1', 'sample2R2'])]
        for (where, exp) in exps:
            dir_fmt = filter_samples(self.sample_paired, self.md_paired_none,
                                     where, True)
            self._assert_paired_contains(dir_fmt, exp)

    def test_filter_paired_id_not_present(self):
        wheres = [None, "Study='A'", "Study='A' OR Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, '1'):
                filter_samples(self.sample_paired, self.md_paired_none, where,
                               False)

    def test_filter_paired_none_no_filter(self):
        wheres = ["Study='B'", "Study='A' AND Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, 'No'):
                filter_samples(self.sample_paired, self.md_paired_none, where,
                               False)

    def test_filter_paired_none_exclude_no_filter(self):
        wheres = ["Study='B'", "Study='A' AND Study='B'"]
        for where in wheres:
            with self.assertRaisesRegex(ValueError, 'No'):
                filter_samples(self.sample_paired, self.md_paired_none, where,
                               True)


if __name__ == '__main__':
    unittest.main()
