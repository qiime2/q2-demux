# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.per_sample_sequences import FastqGzFormat
import qiime2.plugin.model as model


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
