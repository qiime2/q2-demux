# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import ValidationError

from q2_demux._format import ErrorCorrectionDetailsFmt


class TestErrorCorrectionDetailsFmt(TestPluginBase):
    package = 'q2_demux.tests'

    def test_validate_positive(self):
        fp = self.get_data_path('error_correction_details/positive.tsv')
        # Should just work
        ErrorCorrectionDetailsFmt(fp, mode='r').validate()

    def test_validate_invalid_format(self):
        fp = self.get_data_path('error_correction_details/invalid.tsv')
        with self.assertRaisesRegex(ValidationError,
                                    'Failed to locate header.'):
            ErrorCorrectionDetailsFmt(fp, mode='r').validate()

    def test_validate_missing_columns(self):
        fp = self.get_data_path('error_correction_details/missing_columns.tsv')
        with self.assertRaisesRegex(ValidationError,
                                    'barcode-corrected.*is not a column'):
            ErrorCorrectionDetailsFmt(fp, mode='r').validate()
