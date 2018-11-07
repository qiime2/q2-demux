# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import skbio
import numpy.testing as npt
import numpy as np

import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.per_sample_sequences import (
    FastqGzFormat, CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt)
from q2_demux._demux import (BarcodeSequenceFastqIterator,
                             BarcodePairedSequenceFastqIterator)
from q2_demux import (emp_single, emp_paired,
                      subsample_single, subsample_paired)


class SubsampleTests(TestPluginBase):
    # this functionality is derived from test_demux.EmpTestingUtils
    package = 'q2_demux.tests'

    def _decode_qual_to_phred(self, qual_str):
        # this function is adapted from scikit-bio
        qual = np.fromstring(qual_str, dtype=np.uint8) - 33
        return qual

    def _validate_fastq_subsampled(self, fastq, sequences, indices):
        # used for tests where only some input seqs must be present in the
        # output (i.e., some subsampling has been performed)
        output_seqs = skbio.io.read(fastq, format='fastq', phred_offset=33,
                                    compression='gzip', constructor=skbio.DNA)
        output_seqs = list(output_seqs)

        input_seqs = set()
        for idx in indices:
            sequence = sequences[idx]
            input_seqs.add((sequence[0], sequence[1], sequence[2],
                           tuple(self._decode_qual_to_phred(sequence[3]))))

        # the number of output sequences is less than or equal to the
        # number of input sequences
        self.assertTrue(len(output_seqs) <= len(input_seqs))

        # confirm that any output seqs are found in the input seqs, and
        # the header, sequence, and quality all perfectly match the input
        for e in output_seqs:
            header_line = '@%s %s' % (e.metadata['id'],
                                      e.metadata['description'])
            output_seq = (header_line,
                          str(e),
                          '+',
                          tuple(e.positional_metadata['quality']))
            self.assertTrue(output_seq in input_seqs)

        # return the number of output sequences, so that can be used in
        # other tests
        return len(output_seqs)


class SubsampleSingleTests(SubsampleTests):

    def setUp(self):
        super().setUp()

        demuxed = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path('subsample_single_end'), mode='r')
        self.demux_data = transform(
            demuxed, to_type=SingleLanePerSampleSingleEndFastqDirFmt)

    def test_no_subsample(self):
        actual = subsample_single(self.demux_data, fraction=1.0)

        # five forward sample files
        forward_fastq = [
            view for path, view in actual.sequences.iter_views(FastqGzFormat)
            if 'R1_001.fastq' in path.name]
        self.assertEqual(len(forward_fastq), 5)

        # FORWARD:
        fwd_subsampled_sequences = 0
        # sequences in sample1 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[0].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[2].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[3].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[4].open(), self.sequences, [6, 8, 9])

        # all input sequences are included in output
        self.assertEqual(fwd_subsampled_sequences, 11)

    def test_subsample(self):
        actual = subsample_single(self.demux_data, fraction=0.5)

        # five forward sample files
        forward_fastq = [
            view for path, view in actual.sequences.iter_views(FastqGzFormat)
            if 'R1_001.fastq' in path.name]
        self.assertEqual(len(forward_fastq), 5)

        # FORWARD:
        fwd_subsampled_sequences = 0
        # sequences in sample1 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[0].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[2].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[3].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[4].open(), self.sequences, [6, 8, 9])
        # some sequences have been removed - this could occasionally fail,
        # but the frequency of that should be ~ 2 * 0.5 ** 11
        self.assertTrue(0 < fwd_subsampled_sequences < 11)

    def test_correct_output_files_on_small_subsample(self):
        # some or all of the output files are likely to be empty, but they
        # should still be present and in the manifest
        actual = subsample_single(self.demux_data, fraction=0.00001)

        # five forward sample files
        forward_fastq = [
            view for path, view in actual.sequences.iter_views(FastqGzFormat)
            if 'R1_001.fastq' in path.name]
        self.assertEqual(len(forward_fastq), 5)


class SubsamplePairedTests(SubsampleTests):

    def setUp(self):
        super().setUp()
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
        bsi = BarcodePairedSequenceFastqIterator(
            self.barcodes, self.forward, self.reverse)

        barcode_map = pd.Series(
            ['AAAA', 'AACC', 'TTAA', 'GGAA', 'CGGC'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        self.demux_data = emp_paired(bsi, barcode_map)

    def test_no_subsample(self):
        actual = subsample_paired(self.demux_data, fraction=1.0)

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
        fwd_subsampled_sequences = 0
        # sequences in sample1 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[0].open(), self.forward, [0, 5])

        # sequences in sample2 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[1].open(), self.forward, [2, 4])

        # sequences in sample3 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[2].open(), self.forward, [1, 3])

        # sequences in sample4 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[3].open(), self.forward, [7, 10])

        # sequences in sample5 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[4].open(), self.forward, [6, 8, 9])

        # all input sequences are included in output
        self.assertEqual(fwd_subsampled_sequences, 11)

        # REVERSE:
        rev_subsampled_sequences = 0
        # sequences in sample1 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[0].open(), self.reverse, [0, 5])

        # sequences in sample2 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[1].open(), self.reverse, [2, 4])

        # sequences in sample3 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[2].open(), self.reverse, [1, 3])

        # sequences in sample4 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[3].open(), self.reverse, [7, 10])

        # sequences in sample5 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[4].open(), self.reverse, [6, 8, 9])

        # all input sequences are included in output
        self.assertEqual(rev_subsampled_sequences, 11)

    def test_subsample(self):
        actual = subsample_paired(self.demux_data, fraction=0.5)

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
        fwd_subsampled_sequences = 0
        # sequences in sample1 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[0].open(), self.forward, [0, 5])

        # sequences in sample2 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[1].open(), self.forward, [2, 4])

        # sequences in sample3 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[2].open(), self.forward, [1, 3])

        # sequences in sample4 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[3].open(), self.forward, [7, 10])

        # sequences in sample5 are correct
        fwd_subsampled_sequences += self._validate_fastq_subsampled(
            forward_fastq[4].open(), self.forward, [6, 8, 9])
        # some sequences have been removed - this could occasionally fail,
        # but the frequency of that should be ~ 2 * 0.5 ** 11
        self.assertTrue(0 < fwd_subsampled_sequences < 11)

        # REVERSE:
        rev_subsampled_sequences = 0
        # sequences in sample1 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[0].open(), self.reverse, [0, 5])

        # sequences in sample2 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[1].open(), self.reverse, [2, 4])

        # sequences in sample3 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[2].open(), self.reverse, [1, 3])

        # sequences in sample4 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[3].open(), self.reverse, [7, 10])

        # sequences in sample5 are correct
        rev_subsampled_sequences += self._validate_fastq_subsampled(
            reverse_fastq[4].open(), self.reverse, [6, 8, 9])
        # some sequences have been removed - this could occasionally fail,
        # but the frequency of that should be ~ 2 * 0.5 ** 11
        self.assertTrue(0 < rev_subsampled_sequences < 11)

        self.assertEqual(fwd_subsampled_sequences, rev_subsampled_sequences)

    def test_correct_output_files_on_small_subsample(self):
        # some or all of the output files are likely to be empty, but they
        # should still be present and in the manifest
        actual = subsample_paired(self.demux_data, fraction=0.00001)

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


if __name__ == '__main__':
    unittest.main()
