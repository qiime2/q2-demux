# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import unittest.mock as mock
import os.path
import tempfile

import pandas as pd
import skbio
import qiime2
import numpy as np
import numpy.testing as npt

from q2_demux._demux import (BarcodeSequenceFastqIterator,
                             BarcodePairedSequenceFastqIterator)
from q2_demux import emp_single, emp_paired, summarize
from q2_types.per_sample_sequences import (
    FastqGzFormat, FastqManifestFormat, YamlFormat)
from q2_demux._summarize._visualizer import _PlotQualView


class BarcodeSequenceFastqIteratorTests(unittest.TestCase):

    def test_valid(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        for i, (barcode, sequence) in enumerate(bsi):
            self.assertEqual(barcode, barcodes[i])
            self.assertEqual(sequence, sequences[i])

    def test_too_few_barcodes(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_too_few_sequences(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_id(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s5/1 abc/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_description(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abd/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_handles_slashes_in_id(self):
        # mismatch is detected as being before the last slash, even if there
        # is more than one slash
        barcodes = [('@s1/2/2 abc/2', 'AAAA', '+', 'YYYY')]
        sequences = [('@s1/1/1 abc/1', 'GGG', '+', 'YYY')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_handles_slashes_in_description(self):
        # mismatch is detected as being before the last slash, even if there
        # is more than one slash
        barcodes = [('@s1/2 a/2/2', 'AAAA', '+', 'YYYY')]
        sequences = [('@s1/1 a/1/1', 'GGG', '+', 'YYY')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_no_description(self):
        barcodes = [('@s1/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1', 'GGG', '+', 'YYY'),
                     ('@s2/1', 'CCC', '+', 'PPP'),
                     ('@s3/1', 'AAA', '+', 'PPP'),
                     ('@s4/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        for i, (barcode, sequence) in enumerate(bsi):
            self.assertEqual(barcode, barcodes[i])
            self.assertEqual(sequence, sequences[i])

    def test_only_one_description(self):
        barcodes = [('@s1/2 abc', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1', 'GGG', '+', 'YYY'),
                     ('@s2/1', 'CCC', '+', 'PPP'),
                     ('@s3/1', 'AAA', '+', 'PPP'),
                     ('@s4/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

        barcodes = [('@s1/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)


class EmpTestingUtils:
    def _decode_qual_to_phred(self, qual_str):
        # this function is adapted from scikit-bio
        qual = np.fromstring(qual_str, dtype=np.uint8) - 33
        return qual

    def _compare_sequence_to_record(self, sequence, fields):
        header_line = ' '.join([sequence.metadata['id'],
                                sequence.metadata['description']])
        self.assertEqual(fields[0][1:], header_line)
        self.assertEqual(fields[1], str(sequence))
        npt.assert_array_equal(self._decode_qual_to_phred(fields[3]),
                               sequence.positional_metadata['quality'])

    def _compare_manifests(self, act_manifest, exp_manifest):
        # strip comment lines before comparing
        act_manifest = [l for l in act_manifest if not l.startswith('#')]
        self.assertEqual(act_manifest, exp_manifest)

    def _validate_sample_fastq(self, fastq, sequences, indices):
        seqs = skbio.io.read(fastq, format='fastq', phred_offset=33,
                             compression='gzip', constructor=skbio.DNA)
        seqs = list(seqs)
        self.assertEqual(len(seqs), len(indices))
        for idx, i in enumerate(indices):
            self._compare_sequence_to_record(seqs[idx], sequences[i])


class EmpSingleTests(unittest.TestCase, EmpTestingUtils):

    def setUp(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'CGGC', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'GGAA', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'CGGC', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'CGGC', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'GGAA', '+', 'PPPP')]

        self.sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                          ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                          ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                          ('@s4/1 abc/1', 'TTT', '+', 'PPP'),
                          ('@s5/1 abc/1', 'ATA', '+', 'PPP'),
                          ('@s6/1 abc/1', 'TAT', '+', 'PPP'),
                          ('@s7/1 abc/1', 'CGC', '+', 'PPP'),
                          ('@s8/1 abc/1', 'GCG', '+', 'PPP'),
                          ('@s9/1 abc/1', 'ACG', '+', 'PPP'),
                          ('@s10/1 abc/1', 'GCA', '+', 'PPP'),
                          ('@s11/1 abc/1', 'TGA', '+', 'PPP')]
        self.bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)

        barcode_map = pd.Series(['AAAA', 'AACC', 'TTAA', 'GGAA', 'CGGC'],
                                index=['sample1', 'sample2', 'sample3',
                                       'sample4', 'sample5'])
        self.barcode_map = qiime2.MetadataCategory(barcode_map)

    def test_valid(self):
        actual = emp_single(self.bsi, self.barcode_map)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    @mock.patch('q2_demux._demux.OPEN_FH_LIMIT', 3)
    def test_valid_small_open_fh_limit(self):
        self.test_valid()

    def test_variable_length_barcodes(self):
        barcodes = pd.Series(['AAA', 'AACC'], index=['sample1', 'sample2'])
        barcodes = qiime2.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp_single(self.bsi, barcodes)

    def test_duplicate_barcodes(self):
        barcodes = pd.Series(['AACC', 'AACC'], index=['sample1', 'sample2'])
        barcodes = qiime2.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp_single(self.bsi, barcodes)

    def test_no_matched_barcodes(self):
        barcodes = pd.Series(['CCCC', 'GGCC'], index=['sample1', 'sample2'])
        barcodes = qiime2.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp_single(self.bsi, barcodes)

    def test_rev_comp_mapping_barcodes(self):
        barcodes = pd.Series(['TTTT', 'GGTT', 'TTAA', 'TTCC', 'GCCG'],
                             index=['sample1', 'sample2', 'sample3', 'sample4',
                                    'sample5'])
        barcodes = qiime2.MetadataCategory(barcodes)
        actual = emp_single(self.bsi, barcodes, rev_comp_mapping_barcodes=True)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_rev_comp_barcodes(self):
        barcodes = [('@s1/2 abc/2', 'TTTT', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'TTTT', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'TTCC', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'TTCC', '+', 'PPPP')]
        bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)
        actual = emp_single(bsi, self.barcode_map, rev_comp_barcodes=True)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_barcode_trimming(self):
        # these barcodes are longer then the ones in the mapping file, so
        # only the first barcode_length bases should be read
        barcodes = [('@s1/2 abc/2', 'AAAAG', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAAG', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACCG', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAAG', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'AACCG', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'AAAAG', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'GGAAG', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'GGAAG', '+', 'PPPP')]
        bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)
        actual = emp_single(bsi, self.barcode_map)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)


class EmpPairedTests(unittest.TestCase, EmpTestingUtils):
    def setUp(self):
        self.barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                         ('@s2/2 abc/2', 'TTAA', '+', 'PPPP'),
                         ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                         ('@s4/2 abc/2', 'TTAA', '+', 'PPPP'),
                         ('@s5/2 abc/2', 'AACC', '+', 'PPPP'),
                         ('@s6/2 abc/2', 'AAAA', '+', 'PPPP'),
                         ('@s7/2 abc/2', 'CGGC', '+', 'PPPP'),
                         ('@s8/2 abc/2', 'GGAA', '+', 'PPPP'),
                         ('@s9/2 abc/2', 'CGGC', '+', 'PPPP'),
                         ('@s10/2 abc/2', 'CGGC', '+', 'PPPP'),
                         ('@s11/2 abc/2', 'GGAA', '+', 'PPPP')]

        self.forward = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                        ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                        ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                        ('@s4/1 abc/1', 'TTT', '+', 'PPP'),
                        ('@s5/1 abc/1', 'ATA', '+', 'PPP'),
                        ('@s6/1 abc/1', 'TAT', '+', 'PPP'),
                        ('@s7/1 abc/1', 'CGC', '+', 'PPP'),
                        ('@s8/1 abc/1', 'GCG', '+', 'PPP'),
                        ('@s9/1 abc/1', 'ACG', '+', 'PPP'),
                        ('@s10/1 abc/1', 'GCA', '+', 'PPP'),
                        ('@s11/1 abc/1', 'TGA', '+', 'PPP')]

        self.reverse = [('@s1/1 abc/1', 'CCC', '+', 'YYY'),
                        ('@s2/1 abc/1', 'GGG', '+', 'PPP'),
                        ('@s3/1 abc/1', 'TTT', '+', 'PPP'),
                        ('@s4/1 abc/1', 'AAA', '+', 'PPP'),
                        ('@s5/1 abc/1', 'TAT', '+', 'PPP'),
                        ('@s6/1 abc/1', 'ATA', '+', 'PPP'),
                        ('@s7/1 abc/1', 'GCG', '+', 'PPP'),
                        ('@s8/1 abc/1', 'CGC', '+', 'PPP'),
                        ('@s9/1 abc/1', 'CGT', '+', 'PPP'),
                        ('@s10/1 abc/1', 'TGC', '+', 'PPP'),
                        ('@s11/1 abc/1', 'TCA', '+', 'PPP')]

        self.bpsi = BarcodePairedSequenceFastqIterator(
            self.barcodes, self.forward, self.reverse)

        barcode_map = pd.Series(['AAAA', 'AACC', 'TTAA', 'GGAA', 'CGGC'],
                                index=['sample1', 'sample2', 'sample3',
                                       'sample4', 'sample5'])
        self.barcode_map = qiime2.MetadataCategory(barcode_map)

    def check_valid(self, *args, **kwargs):
        actual = emp_paired(*args, **kwargs)

        # five forward sample files
        forward_fastq = [
            view for path, view in actual.sequences.iter_views(FastqGzFormat)
            if 'R1_001.fastq' in path.name]
        self.assertEqual(len(forward_fastq), 5)

        # five reverse sample files
        reverse_fastq = [
            view for path, view in actual.sequences.iter_views(FastqGzFormat)
            if 'R2_001.fastq' in path.name]
        self.assertEqual(len(reverse_fastq), 5)

        # FORWARD:
        # sequences in sample1 are correct
        self._validate_sample_fastq(
            forward_fastq[0].open(), self.forward, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            forward_fastq[1].open(), self.forward, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            forward_fastq[2].open(), self.forward, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            forward_fastq[3].open(), self.forward, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            forward_fastq[4].open(), self.forward, [6, 8, 9])

        # REVERSE:
        # sequences in sample1 are correct
        self._validate_sample_fastq(
            reverse_fastq[0].open(), self.reverse, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            reverse_fastq[1].open(), self.reverse, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            reverse_fastq[2].open(), self.reverse, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            reverse_fastq[3].open(), self.reverse, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            reverse_fastq[4].open(), self.reverse, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample1,sample1_1_L001_R2_001.fastq.gz,reverse\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R2_001.fastq.gz,reverse\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R2_001.fastq.gz,reverse\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R2_001.fastq.gz,reverse\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R2_001.fastq.gz,reverse\n']

        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_valid(self):
        self.check_valid(self.bpsi, self.barcode_map)

    @mock.patch('q2_demux._demux.OPEN_FH_LIMIT', 6)
    def test_valid_small_open_fh_limit(self):
        self.test_valid()

    def test_variable_length_barcodes(self):
        barcodes = pd.Series(['AAA', 'AACC'], index=['sample1', 'sample2'])
        barcodes = qiime2.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp_paired(self.bpsi, barcodes)

    def test_duplicate_barcodes(self):
        barcodes = pd.Series(['AACC', 'AACC'], index=['sample1', 'sample2'])
        barcodes = qiime2.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp_paired(self.bpsi, barcodes)

    def test_no_matched_barcodes(self):
        barcodes = pd.Series(['CCCC', 'GGCC'], index=['sample1', 'sample2'])
        barcodes = qiime2.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp_paired(self.bpsi, barcodes)

    def test_rev_comp_mapping_barcodes(self):
        barcodes = pd.Series(['TTTT', 'GGTT', 'TTAA', 'TTCC', 'GCCG'],
                             index=['sample1', 'sample2', 'sample3', 'sample4',
                                    'sample5'])
        barcodes = qiime2.MetadataCategory(barcodes)
        self.check_valid(self.bpsi, barcodes, rev_comp_mapping_barcodes=True)

    def test_rev_comp_barcodes(self):
        barcodes = [('@s1/2 abc/2', 'TTTT', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'TTTT', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'TTCC', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'TTCC', '+', 'PPPP')]
        bpsi = BarcodePairedSequenceFastqIterator(
            barcodes, self.forward, self.reverse)
        self.check_valid(bpsi, self.barcode_map, rev_comp_barcodes=True)

    def test_barcode_trimming(self):
        # these barcodes are longer then the ones in the mapping file, so
        # only the first barcode_length bases should be read
        barcodes = [('@s1/2 abc/2', 'AAAAG', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAAG', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACCG', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAAG', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'AACCG', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'AAAAG', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'GGAAG', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'GGAAG', '+', 'PPPP')]
        bpsi = BarcodePairedSequenceFastqIterator(
            barcodes, self.forward, self.reverse)
        self.check_valid(bpsi, self.barcode_map)


class SummarizeTests(unittest.TestCase):

    def setUp(self):
        self.barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                         ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                         ('@s3/2 abc/2', 'AAAA', '+', 'PPPP'),
                         ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        self.sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                          ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                          ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                          ('@s4/1 abc/1', 'TTT', '+', 'PPP')]

    def test_basic(self):
        bsi = BarcodeSequenceFastqIterator(self.barcodes, self.sequences)

        barcode_map = pd.Series(['AAAA', 'AACC'], index=['sample1', 'sample2'])
        barcode_map = qiime2.MetadataCategory(barcode_map)

        demux_data = emp_single(bsi, barcode_map)
        # test that an index.html file is created and that it has size > 0
        with tempfile.TemporaryDirectory() as output_dir:
            # TODO: Remove _PlotQualView wrapper
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=2)
            self.assertTrue(result is None)
            index_fp = os.path.join(output_dir, 'overview.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.getsize(index_fp) > 0)
            csv_fp = os.path.join(output_dir, 'per-sample-fastq-counts.csv')
            self.assertTrue(os.path.exists(csv_fp))
            self.assertTrue(os.path.getsize(csv_fp) > 0)
            pdf_fp = os.path.join(output_dir, 'demultiplex-summary.pdf')
            self.assertTrue(os.path.exists(pdf_fp))
            self.assertTrue(os.path.getsize(pdf_fp) > 0)
            png_fp = os.path.join(output_dir, 'demultiplex-summary.png')
            self.assertTrue(os.path.exists(png_fp))
            self.assertTrue(os.path.getsize(png_fp) > 0)
            with open(index_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<td>Minimum:</td><td>1</td>', html)
                self.assertIn('<td>Maximum:</td><td>3</td>', html)

    def test_single_sample(self):
        bsi = BarcodeSequenceFastqIterator(self.barcodes[:1],
                                           self.sequences[:1])

        barcode_map = pd.Series(['AAAA'], index=['sample1'])
        barcode_map = qiime2.MetadataCategory(barcode_map)

        demux_data = emp_single(bsi, barcode_map)
        # test that an index.html file is created and that it has size > 0
        with tempfile.TemporaryDirectory() as output_dir:
            # TODO: Remove _PlotQualView wrapper
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=1)
            self.assertTrue(result is None)
            index_fp = os.path.join(output_dir, 'overview.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.getsize(index_fp) > 0)
            csv_fp = os.path.join(output_dir, 'per-sample-fastq-counts.csv')
            self.assertTrue(os.path.exists(csv_fp))
            self.assertTrue(os.path.getsize(csv_fp) > 0)
            pdf_fp = os.path.join(output_dir, 'demultiplex-summary.pdf')
            self.assertFalse(os.path.exists(pdf_fp))
            png_fp = os.path.join(output_dir, 'demultiplex-summary.png')
            self.assertFalse(os.path.exists(png_fp))
            with open(index_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<td>Minimum:</td><td>1</td>', html)
                self.assertIn('<td>Maximum:</td><td>1</td>', html)

    def test_paired_end(self):
        barcodes = self.barcodes[:3]

        forward = self.sequences[:3]

        reverse = [('@s1/1 abc/1', 'CCC', '+', 'YYY'),
                   ('@s2/1 abc/1', 'GGG', '+', 'PPP'),
                   ('@s3/1 abc/1', 'TTT', '+', 'PPP')]

        bpsi = BarcodePairedSequenceFastqIterator(barcodes, forward, reverse)

        barcode_map = pd.Series(['AAAA', 'AACC', 'TTAA'],
                                index=['sample1', 'sample2', 'sample3'])
        barcode_map = qiime2.MetadataCategory(barcode_map)

        demux_data = emp_paired(bpsi, barcode_map)
        with tempfile.TemporaryDirectory() as output_dir:
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=True), n=2)
            self.assertTrue(result is None)
            plot_fp = os.path.join(output_dir, 'quality-plot.html')
            with open(plot_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<h5 class="text-center">Forward Reads</h5>',
                              html)
                self.assertIn('<h5 class="text-center">Reverse Reads</h5>',
                              html)

    def test_subsample_higher_than_seqs_count(self):
        barcodes = self.barcodes[:1]

        sequences = self.sequences[:1]
        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)

        barcode_map = pd.Series(['AAAA'], index=['sample1'])
        barcode_map = qiime2.MetadataCategory(barcode_map)

        demux_data = emp_single(bsi, barcode_map)
        with tempfile.TemporaryDirectory() as output_dir:
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=50)
            self.assertTrue(result is None)
            plot_fp = os.path.join(output_dir, 'quality-plot.html')
            with open(plot_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<strong>Warning:</strong>', html)

    def test_phred_score_out_of_range(self):
        barcodes = self.barcodes[:3]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'jjj'),
                     ('@s2/1 abc/1', 'CCC', '+', 'iii'),
                     ('@s3/1 abc/1', 'AAA', '+', 'hhh')]
        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)

        barcode_map = pd.Series(['AAAA', 'AACC', 'TTAA'],
                                index=['sample1', 'sample2', 'sample3'])
        barcode_map = qiime2.MetadataCategory(barcode_map)

        demux_data = emp_single(bsi, barcode_map)
        with tempfile.TemporaryDirectory() as output_dir:
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=50)
            self.assertTrue(result is None)
            plot_fp = os.path.join(output_dir, 'quality-plot.html')
            with open(plot_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<strong>Danger:</strong>', html)

    def test_inconsistent_sequence_length_single(self):
        sequences = [('@s1/1 abc/1', 'GGGGGGG', '+', 'YYYYYYY'),
                     ('@s2/1 abc/1', 'CCCCC', '+', 'PPPPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc/1', 'T', '+', 'P')]
        bsi = BarcodeSequenceFastqIterator(self.barcodes, sequences)

        barcode_map = pd.Series(['AAAA', 'AACC'], index=['sample1', 'sample2'])
        barcode_map = qiime2.MetadataCategory(barcode_map)

        demux_data = emp_single(bsi, barcode_map)
        with tempfile.TemporaryDirectory() as output_dir:
            # TODO: Remove _PlotQualView wrapper
            summarize(output_dir, _PlotQualView(demux_data,
                                                paired=False), n=4)
            plot_fp = os.path.join(output_dir, 'data.jsonp')
            with open(plot_fp, 'r') as fh:
                jsonp = fh.read()
                self.assertIn('"seqCount":4', jsonp)
                self.assertIn('"minSeqLen":1', jsonp)
                self.assertIn('"n":4', jsonp)

    def test_inconsistent_sequence_length_paired(self):
        forward = [('@s1/1 abc/1', 'G', '+', 'Y'),
                   ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                   ('@s3/1 abc/1', 'AAAAA', '+', 'PPPPP'),
                   ('@s4/1 abc/1', 'TTTTTTT', '+', 'PPPPPPP')]
        reverse = [('@s1/1 abc/1', 'AAAAAAA', '+', 'YYYYYYY'),
                   ('@s2/1 abc/1', 'TTTTT', '+', 'PPPPP'),
                   ('@s3/1 abc/1', 'GGG', '+', 'PPP'),
                   ('@s4/1 abc/1', 'C', '+', 'P')]
        bpsi = BarcodePairedSequenceFastqIterator(self.barcodes, forward,
                                                  reverse)

        barcode_map = pd.Series(['AAAA', 'AACC'], index=['sample1', 'sample2'])
        barcode_map = qiime2.MetadataCategory(barcode_map)

        demux_data = emp_paired(bpsi, barcode_map)
        with tempfile.TemporaryDirectory() as output_dir:
            # TODO: Remove _PlotQualView wrapper
            summarize(output_dir, _PlotQualView(demux_data,
                                                paired=True), n=4)
            plot_fp = os.path.join(output_dir, 'data.jsonp')
            with open(plot_fp, 'r') as fh:
                jsonp = fh.read()
                self.assertIn('"seqCount":4', jsonp)
                self.assertIn('"minSeqLen":1', jsonp)
                self.assertIn('"n":4', jsonp)


if __name__ == '__main__':
    unittest.main()
