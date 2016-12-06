import unittest
import os.path
import tempfile

import pandas as pd
import skbio
import qiime
import numpy as np
import numpy.testing as npt

from q2_demux._demux import BarcodeSequenceFastqIterator
from q2_demux import emp, summarize
from q2_types.per_sample_sequences import (
    FastqGzFormat, FastqManifestFormat, YamlFormat)


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


class EmpTests(unittest.TestCase):

    def setUp(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        self.sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                          ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                          ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                          ('@s4/1 abc/1', 'TTT', '+', 'PPP')]
        self.bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)

        barcode_map = pd.Series(['AAAA', 'AACC'], index=['sample1', 'sample2'])
        self.barcode_map = qiime.MetadataCategory(barcode_map)

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

    def test_valid(self):
        actual = emp(self.bsi, self.barcode_map)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # two per-sample files were written
        self.assertEqual(len(output_fastq), 2)

        # sequences in sample1 are correct
        sample1_seqs = skbio.io.read(output_fastq[0][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample1_seqs = list(sample1_seqs)
        self.assertEqual(len(sample1_seqs), 2)
        self._compare_sequence_to_record(sample1_seqs[0], self.sequences[0])
        self._compare_sequence_to_record(sample1_seqs[1], self.sequences[1])

        # sequences in sample2 are correct
        sample2_seqs = skbio.io.read(output_fastq[1][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample2_seqs = list(sample2_seqs)
        self.assertEqual(len(sample2_seqs), 2)
        self._compare_sequence_to_record(sample2_seqs[0], self.sequences[2])
        self._compare_sequence_to_record(sample2_seqs[1], self.sequences[3])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_2_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_variable_length_barcodes(self):
        barcodes = pd.Series(['AAA', 'AACC'], index=['sample1', 'sample2'])
        barcodes = qiime.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp(self.bsi, barcodes)

    def test_duplicate_barcodes(self):
        barcodes = pd.Series(['AACC', 'AACC'], index=['sample1', 'sample2'])
        barcodes = qiime.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp(self.bsi, barcodes)

    def test_no_matched_barcodes(self):
        barcodes = pd.Series(['CCCC', 'GGCC'], index=['sample1', 'sample2'])
        barcodes = qiime.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp(self.bsi, barcodes)

    def test_rev_comp_mapping_barcodes(self):
        barcodes = pd.Series(['TTTT', 'GGTT'], index=['sample1', 'sample2'])
        barcodes = qiime.MetadataCategory(barcodes)
        actual = emp(self.bsi, barcodes, rev_comp_mapping_barcodes=True)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # two per-sample files were written
        self.assertEqual(len(output_fastq), 2)

        # sequences in sample1 are correct
        sample1_seqs = skbio.io.read(output_fastq[0][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample1_seqs = list(sample1_seqs)
        self.assertEqual(len(sample1_seqs), 2)
        self._compare_sequence_to_record(sample1_seqs[0], self.sequences[0])
        self._compare_sequence_to_record(sample1_seqs[1], self.sequences[1])

        # sequences in sample2 are correct
        sample2_seqs = skbio.io.read(output_fastq[1][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample2_seqs = list(sample2_seqs)
        self.assertEqual(len(sample2_seqs), 2)
        self._compare_sequence_to_record(sample2_seqs[0], self.sequences[2])
        self._compare_sequence_to_record(sample2_seqs[1], self.sequences[3])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_2_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_rev_comp_barcodes(self):
        barcodes = [('@s1/2 abc/2', 'TTTT', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTTT', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'GGTT', '+', 'PPPP')]
        bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)
        actual = emp(bsi, self.barcode_map, rev_comp_barcodes=True)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # two per-sample files were written
        self.assertEqual(len(output_fastq), 2)

        # sequences in sample1 are correct
        sample1_seqs = skbio.io.read(output_fastq[0][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample1_seqs = list(sample1_seqs)
        self.assertEqual(len(sample1_seqs), 2)
        self._compare_sequence_to_record(sample1_seqs[0], self.sequences[0])
        self._compare_sequence_to_record(sample1_seqs[1], self.sequences[1])

        # sequences in sample2 are correct
        sample2_seqs = skbio.io.read(output_fastq[1][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample2_seqs = list(sample2_seqs)
        self.assertEqual(len(sample2_seqs), 2)
        self._compare_sequence_to_record(sample2_seqs[0], self.sequences[2])
        self._compare_sequence_to_record(sample2_seqs[1], self.sequences[3])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_2_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_barcode_trimming(self):
        # these barcodes are longer then the ones in the mapping file, so
        # only the first barcode_length bases should be read
        barcodes = [('@s1/2 abc/2', 'AAAAG', '+', 'YYYYY'),
                    ('@s2/2 abc/2', 'AAAAG', '+', 'PPPPP'),
                    ('@s3/2 abc/2', 'AACCG', '+', 'PPPPP'),
                    ('@s4/2 abc/2', 'AACCG', '+', 'PPPPP')]
        bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)
        actual = emp(bsi, self.barcode_map)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # two per-sample files were written
        self.assertEqual(len(output_fastq), 2)

        # sequences in sample1 are correct
        sample1_seqs = skbio.io.read(output_fastq[0][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample1_seqs = list(sample1_seqs)
        self.assertEqual(len(sample1_seqs), 2)
        self._compare_sequence_to_record(sample1_seqs[0], self.sequences[0])
        self._compare_sequence_to_record(sample1_seqs[1], self.sequences[1])

        # sequences in sample2 are correct
        sample2_seqs = skbio.io.read(output_fastq[1][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample2_seqs = list(sample2_seqs)
        self.assertEqual(len(sample2_seqs), 2)
        self._compare_sequence_to_record(sample2_seqs[0], self.sequences[2])
        self._compare_sequence_to_record(sample2_seqs[1], self.sequences[3])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_2_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)


class SummarizeTests(unittest.TestCase):

    def test_basic(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc/1', 'TTT', '+', 'PPP')]
        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)

        barcode_map = pd.Series(['AAAA', 'AACC'], index=['sample1', 'sample2'])
        barcode_map = qiime.MetadataCategory(barcode_map)

        demux_data = emp(bsi, barcode_map)
        # test that an index.html file is created and that it has size > 0
        with tempfile.TemporaryDirectory() as output_dir:
            result = summarize(output_dir, demux_data)
            self.assertTrue(result is None)
            index_fp = os.path.join(output_dir, 'index.html')
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
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY')]
        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)

        barcode_map = pd.Series(['AAAA'], index=['sample1'])
        barcode_map = qiime.MetadataCategory(barcode_map)

        demux_data = emp(bsi, barcode_map)
        # test that an index.html file is created and that it has size > 0
        with tempfile.TemporaryDirectory() as output_dir:
            result = summarize(output_dir, demux_data)
            self.assertTrue(result is None)
            index_fp = os.path.join(output_dir, 'index.html')
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
