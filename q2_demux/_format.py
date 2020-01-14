# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.per_sample_sequences import FastqGzFormat
import qiime2.plugin.model as model
from qiime2.plugin import ValidationError
import qiime2


# TODO: deprecate this and alias it
class EMPMultiplexedDirFmt(model.DirectoryFormat):
    sequences = model.File(
        r'sequences.fastq.gz', format=FastqGzFormat)

    barcodes = model.File(
        r'barcodes.fastq.gz', format=FastqGzFormat)


# The new cannonical name for EMPMultiplexedDirFmt
class EMPSingleEndDirFmt(EMPMultiplexedDirFmt):
    pass  # contents inherited


class EMPPairedEndDirFmt(model.DirectoryFormat):
    forward = model.File(
        r'forward.fastq.gz', format=FastqGzFormat)

    reverse = model.File(
        r'reverse.fastq.gz', format=FastqGzFormat)

    barcodes = model.File(
        r'barcodes.fastq.gz', format=FastqGzFormat)


# Originally called EMPMultiplexedSingleEndDirFmt, rename was possible as no
# artifacts where created with this view, it is just for import.
class EMPSingleEndCasavaDirFmt(model.DirectoryFormat):
    # TODO: generalize this with a regex when we have validation in place for
    # model.FileCollections. The file names are currently designed more
    # specificially for handling MiSeq data.
    sequences = model.File(
        r'Undetermined_S0_L001_R1_001.fastq.gz', format=FastqGzFormat)

    barcodes = model.File(
        r'Undetermined_S0_L001_I1_001.fastq.gz', format=FastqGzFormat)


class EMPPairedEndCasavaDirFmt(model.DirectoryFormat):
    forward = model.File(
        r'Undetermined_S0_L001_R1_001.fastq.gz', format=FastqGzFormat)

    reverse = model.File(
        r'Undetermined_S0_L001_R2_001.fastq.gz', format=FastqGzFormat)

    barcodes = model.File(
        r'Undetermined_S0_L001_I1_001.fastq.gz', format=FastqGzFormat)


class ErrorCorrectionDetailsFmt(model.TextFileFormat):
    METADATA_COLUMNS = {
        'sample',
        'barcode-sequence-id',
        'barcode-uncorrected',
        'barcode-corrected',
        'barcode-errors',
    }

    def _validate_(self, level):
        try:
            md = qiime2.Metadata.load(str(self))
        except qiime2.metadata.MetadataFileError as md_exc:
            raise ValidationError(md_exc) from md_exc

        for column in sorted(self.METADATA_COLUMNS):
            try:
                md.get_column(column)
            except ValueError as md_exc:
                raise ValidationError(md_exc) from md_exc


ErrorCorrectionDetailsDirFmt = model.SingleFileDirectoryFormat(
    'ErrorCorrectionDetailsDirFmt', 'details.tsv', ErrorCorrectionDetailsFmt)
