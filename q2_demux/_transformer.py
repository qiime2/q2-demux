# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import shutil

from q2_types.per_sample_sequences import FastqGzFormat

from .plugin_setup import plugin
from ._demux import (BarcodeSequenceFastqIterator,
                     BarcodePairedSequenceFastqIterator, _read_fastq_seqs)
from ._format import (
    EMPMultiplexedDirFmt, EMPMultiplexedSingleEndDirFmt,
    PairedEMPMultiplexedDirFmt, EMPMultiplexedPairedEndDirFmt)


@plugin.register_transformer
def _1(dirfmt: EMPMultiplexedDirFmt) -> BarcodeSequenceFastqIterator:
    barcode_generator = _read_fastq_seqs(
        str(dirfmt.barcodes.view(FastqGzFormat)))
    sequence_generator = _read_fastq_seqs(
        str(dirfmt.sequences.view(FastqGzFormat)))
    result = BarcodeSequenceFastqIterator(barcode_generator,
                                          sequence_generator)
    # ensure that dirfmt stays in scope as long as result does so these
    # generators will work.
    result.__dirfmt = dirfmt
    return result


@plugin.register_transformer
def _2(dirfmt: EMPMultiplexedSingleEndDirFmt) -> EMPMultiplexedDirFmt:
    # TODO: revisit this API to simpify defining transformers
    result = EMPMultiplexedDirFmt().path

    sequences_fp = str(result / 'sequences.fastq.gz')
    barcodes_fp = str(result / 'barcodes.fastq.gz')
    shutil.copyfile(str(dirfmt.sequences.view(FastqGzFormat)), sequences_fp)
    shutil.copyfile(str(dirfmt.barcodes.view(FastqGzFormat)), barcodes_fp)

    return result


@plugin.register_transformer
def _3(dirfmt: EMPMultiplexedPairedEndDirFmt) -> PairedEMPMultiplexedDirFmt:
    result = EMPMultiplexedDirFmt()
    root = result.path

    forward_fp = str(root / 'forward.fastq.gz')
    reverse_fp = str(root / 'reverse.fastq.gz')
    barcodes_fp = str(root / 'barcodes.fastq.gz')
    shutil.copyfile(str(dirfmt.forward.view(FastqGzFormat)), forward_fp)
    shutil.copyfile(str(dirfmt.reverse.view(FastqGzFormat)), reverse_fp)
    shutil.copyfile(str(dirfmt.barcodes.view(FastqGzFormat)), barcodes_fp)

    return result


@plugin.register_transformer
def _5(dirfmt: PairedEMPMultiplexedDirFmt) \
        -> BarcodePairedSequenceFastqIterator:

    barcode_generator = _read_fastq_seqs(
        str(dirfmt.barcodes.view(FastqGzFormat)))
    forward_generator = _read_fastq_seqs(
        str(dirfmt.forward.view(FastqGzFormat)))
    reverse_generator = _read_fastq_seqs(
        str(dirfmt.reverse.view(FastqGzFormat)))
    result = BarcodePairedSequenceFastqIterator(barcode_generator,
                                                forward_generator,
                                                reverse_generator)
    # ensure that dirfmt stays in scope as long as result does so these
    # generators will work.
    result.__dirfmt = dirfmt
    return result
