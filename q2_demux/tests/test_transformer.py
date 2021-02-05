# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile

import pandas as pd
import pandas.testing as pdt

from q2_demux._format import (EMPSingleEndDirFmt,
                              EMPSingleEndCasavaDirFmt,
                              ErrorCorrectionDetailsFmt)
from q2_demux._demux import BarcodeSequenceFastqIterator
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import ValidationError
import qiime2


class TestTransformers(TestPluginBase):
    package = 'q2_demux.tests'

    def setUp(self):
        # TODO generalize plugin lookup when ported to framework. This code
        # is adapted from the base class.
        try:
            from q2_demux.plugin_setup import plugin
        except ImportError:
            self.fail("Could not import plugin object.")

        self.plugin = plugin

        # TODO use qiime temp dir when ported to framework, and when the
        # configurable temp dir exists
        self.temp_dir = tempfile.TemporaryDirectory(
            prefix='q2-demux-test-temp-')

    def test_emp_multiplexed_format_barcode_sequence_iterator(self):
        transformer = self.get_transformer(EMPSingleEndDirFmt,
                                           BarcodeSequenceFastqIterator)
        dirname = 'emp_multiplexed'
        dirpath = self.get_data_path(dirname)
        bsi = transformer(EMPSingleEndDirFmt(dirpath, mode='r'))
        bsi = list(bsi)
        self.assertEqual(len(bsi), 250)
        self.assertEqual(
            bsi[0][0],
            ('@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0',
             'TTAGGCATCTCG',
             '+',
             'B@@FFFFFHHHH'))
        self.assertEqual(
            bsi[0][1],
            ('@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0',
             'GCTTAGGGATTTTATTGTTATCAGGGTTAATCGTGCCAAGAAAAGCGGCATGGTCAATATAAC'
             'CAGTAGTGTTAACAGTCGGGAGAGGAGTGGCATTAACACCATCCTTCATGAACTTAATCCACT'
             'GTTCACCATAAACGTGACGATGAGG',
             '+',
             'C@CFFFFFHHFHHGIJJ?FFHEIIIIHGEIIFHGIIJHGIGBGB?DHIIJJJJCFCHIEGIGG'
             'HGFAEDCEDBCCEEA.;>?BB=288A?AB709@:3:A:C88CCD@CC444@>>34>>ACC:?C'
             'CD<CDCA>A@A>:<?B@?<((2(>?'))

    def test_emp_se_multiplexed_format_barcode_sequence_iterator(self):
        transformer1 = self.get_transformer(EMPSingleEndCasavaDirFmt,
                                            EMPSingleEndDirFmt)
        transformer2 = self.get_transformer(EMPSingleEndDirFmt,
                                            BarcodeSequenceFastqIterator)
        dirname = 'emp_multiplexed_single_end'
        dirpath = self.get_data_path(dirname)
        emp_demultiplexed = \
            transformer1(EMPSingleEndCasavaDirFmt(dirpath, mode='r'))
        bsi = transformer2(EMPSingleEndDirFmt(emp_demultiplexed, mode='r'))
        bsi = list(bsi)
        self.assertEqual(len(bsi), 250)
        self.assertEqual(
            bsi[0][0],
            ('@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0',
             'TTAGGCATCTCG',
             '+',
             'B@@FFFFFHHHH'))
        self.assertEqual(
            bsi[0][1],
            ('@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0',
             'GCTTAGGGATTTTATTGTTATCAGGGTTAATCGTGCCAAGAAAAGCGGCATGGTCAATATAAC'
             'CAGTAGTGTTAACAGTCGGGAGAGGAGTGGCATTAACACCATCCTTCATGAACTTAATCCACT'
             'GTTCACCATAAACGTGACGATGAGG',
             '+',
             'C@CFFFFFHHFHHGIJJ?FFHEIIIIHGEIIFHGIIJHGIGBGB?DHIIJJJJCFCHIEGIGG'
             'HGFAEDCEDBCCEEA.;>?BB=288A?AB709@:3:A:C88CCD@CC444@>>34>>ACC:?C'
             'CD<CDCA>A@A>:<?B@?<((2(>?'))

    def test_invalid(self):
        dirname = 'bad'
        dirpath = self.get_data_path(dirname)
        transformer = self.get_transformer(EMPSingleEndDirFmt,
                                           BarcodeSequenceFastqIterator)
        with self.assertRaises(ValidationError):
            transformer(EMPSingleEndDirFmt(dirpath, mode='r'))

        transformer = self.get_transformer(EMPSingleEndCasavaDirFmt,
                                           EMPSingleEndDirFmt)
        with self.assertRaises(ValidationError):
            transformer(EMPSingleEndCasavaDirFmt(dirpath, 'r'))


class TestErrorCorrectionDetailsFmtTransformers(TestPluginBase):
    package = 'q2_demux.tests'

    def setUp(self):
        super().setUp()

        self.df = pd.DataFrame([
                ['s1', 'seq-1',  'AAC', 'AAA', 2.],
                ['s1', 'seq-4',  'ACA', 'AAA', 20.],
                ['s2', 'seq-5',  'CCA', 'CCC', 1.],
                ['s3', 'seq-50', 'GGT', 'GGG', 1.],
            ],
            columns=['sample', 'barcode-sequence-id', 'barcode-uncorrected',
                     'barcode-corrected', 'barcode-errors'],
            index=pd.Index(['record-1', 'record-2', 'record-3', 'record-4'],
                           name='id'))

        self.serialized = (
            'id\tsample\tbarcode-sequence-id\tbarcode-uncorrected\t'
            'barcode-corrected\tbarcode-errors\n'
            '#q2:types\tcategorical\tcategorical\tcategorical\tcategorical\t'
            'numeric\n'
            'record-1\ts1\tseq-1\tAAC\tAAA\t2\n'
            'record-2\ts1\tseq-4\tACA\tAAA\t20\n'
            'record-3\ts2\tseq-5\tCCA\tCCC\t1\n'
            'record-4\ts3\tseq-50\tGGT\tGGG\t1\n'
        )

    def test_df_to_error_correction_details_fmt(self):
        transformer = self.get_transformer(
            pd.DataFrame, ErrorCorrectionDetailsFmt)
        obs = transformer(self.df)

        with obs.open() as obs:
            self.assertEqual(obs.read(), self.serialized)

    def test_error_correction_details_fmt_to_df(self):
        transformer = self.get_transformer(
            ErrorCorrectionDetailsFmt, pd.DataFrame)
        ff = ErrorCorrectionDetailsFmt()
        with ff.open() as fh:
            fh.write(self.serialized)
        obs = transformer(ff)

        pdt.assert_frame_equal(obs, self.df)

    def test_error_correction_details_fmt_to_metadata(self):
        transformer = self.get_transformer(
            ErrorCorrectionDetailsFmt, qiime2.Metadata)
        ff = ErrorCorrectionDetailsFmt()
        with ff.open() as fh:
            fh.write(self.serialized)
        obs = transformer(ff)

        self.assertEqual(obs, qiime2.Metadata(self.df))


if __name__ == "__main__":
    unittest.main()
