# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile

from q2_demux._format import (EMPSingleEndDirFmt,
                              EMPSingleEndCasavaDirFmt)
from q2_demux._demux import BarcodeSequenceFastqIterator
from qiime2.plugin.testing import TestPluginBase


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
        with self.assertRaises(ValueError):
            transformer(EMPSingleEndDirFmt(dirpath, mode='r'))

        transformer = self.get_transformer(EMPSingleEndCasavaDirFmt,
                                           EMPSingleEndDirFmt)
        with self.assertRaises(ValueError):
            transformer(EMPSingleEndCasavaDirFmt(dirpath, 'r'))


if __name__ == "__main__":
    unittest.main()
