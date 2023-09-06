import pandas as pd

import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_demux import tabulate_counts


class TabulateTests(TestPluginBase):
    package = 'q2_demux.tests'

    def setUp(self):
        super().setUp()

        demuxed_se_1 = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path('tabulate_counts_single_end_1'), mode='r')
        self.demux_se_data_1 = transform(
            demuxed_se_1, to_type=SingleLanePerSampleSingleEndFastqDirFmt)

        demuxed_se_2 = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path('tabulate_counts_single_end_2'), mode='r')
        self.demux_se_data_2 = transform(
            demuxed_se_2, to_type=SingleLanePerSampleSingleEndFastqDirFmt)

        demuxed_pe_1 = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path('tabulate_counts_paired_end_1'), mode='r')
        self.demux_pe_data_1 = transform(
            demuxed_pe_1, to_type=SingleLanePerSamplePairedEndFastqDirFmt)

    def test_tabulate_counts_se(self):
        actual = tabulate_counts([self.demux_se_data_1])

        expected = {'sample1': 2,
                    'sample2': 2,
                    'sample3': 2,
                    'sample4': 2,
                    'sample5': 3}
        expected = pd.Series(expected)
        expected.name = 'Demultiplexed sequence count'
        expected = expected.to_frame()
        expected.index.name = 'sample-id'
        expected = qiime2.Metadata(expected)

        self.assertEqual(actual, expected)

        actual = tabulate_counts([self.demux_se_data_2])

        expected = {'sample6': 2,
                    'sample7': 2}
        expected = pd.Series(expected)
        expected.name = 'Demultiplexed sequence count'
        expected = expected.to_frame()
        expected.index.name = 'sample-id'
        expected = qiime2.Metadata(expected)

        self.assertEqual(actual, expected)

    def test_tabulate_counts_pe(self):
        actual = tabulate_counts([self.demux_pe_data_1])

        expected = {'sample1': 2}
        expected = pd.Series(expected)
        expected.name = 'Demultiplexed sequence count'
        expected = expected.to_frame()
        expected.index.name = 'sample-id'
        expected = qiime2.Metadata(expected)

        self.assertEqual(actual, expected)

    def test_tabulate_counts_multiple(self):
        actual = tabulate_counts([self.demux_se_data_1, self.demux_se_data_2])

        expected = {'sample1': 2,
                    'sample2': 2,
                    'sample3': 2,
                    'sample4': 2,
                    'sample5': 3,
                    'sample6': 2,
                    'sample7': 2}
        expected = pd.Series(expected)
        expected.name = 'Demultiplexed sequence count'
        expected = expected.to_frame()
        expected.index.name = 'sample-id'
        expected = qiime2.Metadata(expected)

        self.assertEqual(actual, expected)

        actual = tabulate_counts([self.demux_pe_data_1, self.demux_se_data_2])

        expected = {'sample1': 2,
                    'sample6': 2,
                    'sample7': 2}
        expected = pd.Series(expected)
        expected.name = 'Demultiplexed sequence count'
        expected = expected.to_frame()
        expected.index.name = 'sample-id'
        expected = qiime2.Metadata(expected)

        self.assertEqual(actual, expected)

    def test_tabulate_counts_error(self):
        with self.assertRaisesRegex(KeyError, 'duplicated.*sample1'):
            tabulate_counts([self.demux_se_data_1, self.demux_pe_data_1])
