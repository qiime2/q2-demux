# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import yaml
import itertools
import collections
import collections.abc
import random
import resource
import re

import skbio
import psutil

import qiime2
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqManifestFormat, YamlFormat)
from ._ecc import GolayDecoder
from ._format import ErrorCorrectionDetailsFmt


FastqHeader = collections.namedtuple('FastqHeader', ['id', 'description'])


class ECDetails:
    COLUMNS = ['id',
               'sample',
               'barcode-sequence-id',
               'barcode-uncorrected',
               'barcode-corrected',
               'barcode-errors']

    def __init__(self, fmt):
        self._fp = open(str(fmt), 'w')
        self._write_header()

    def write(self, parts):
        self._fp.write('\t'.join([str(part) for part in parts]))
        self._fp.write('\n')

    def _write_header(self):
        self.write(self.COLUMNS)


def _read_fastq_seqs(filepath):
    # This function is adapted from @jairideout's SO post:
    # http://stackoverflow.com/a/39302117/3424666
    fh = gzip.open(filepath, 'rt')
    for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
        yield (seq_header.strip(), seq.strip(), qual_header.strip(),
               qual.strip())


def _trim_id(id):
    return id.rsplit('/', 1)[0]


def _trim_description(desc):
    # The first number of ':' seperated description is the read number
    if ':' in desc:
        desc = desc.split(':', 1)[1]
    return desc.rsplit('/', 1)[0]


def _record_to_fastq_header(record):
    tokens = record[0][1:].split(' ', maxsplit=1)
    if len(tokens) == 1:
        id, = tokens
        description = None
    else:
        id, description = tokens

    return FastqHeader(id=id, description=description)


# This is global so that it can be tested without changing the actual ulimits.
# NOTE: UNIX only
OPEN_FH_LIMIT, _ = resource.getrlimit(resource.RLIMIT_NOFILE)


def _maintain_open_fh_count(per_sample_fastqs, paired=False):
    files_to_open = 1 if not paired else 2
    # NOTE: UNIX only
    if psutil.Process().num_fds() + files_to_open < OPEN_FH_LIMIT:
        return

    # currently open per-sample files
    if not paired:
        open_fhs = [fh for fh in per_sample_fastqs.values()
                    if not fh.closed]
    else:
        open_fhs = [fh for fh in per_sample_fastqs.values()
                    if not fh[0].closed]

    # If the number of open files reaches the allotted resources limit,
    # close around 15% of the open files. 15% was chosen because if you
    # only close a single file it will start to have to do it on every loop
    # and on a 35 file benchmark using a hard coded limit of 10 files,
    # only closing one file added an increased runtime of 160%
    n_to_close = round(0.15 * len(open_fhs))
    if paired:
        n_to_close //= 2
    # Never close more than files than are open, also if closing,
    #  close at least the number of files that will need to be opened.
    n_to_close = min(len(open_fhs), max(n_to_close, files_to_open))
    for rand_fh in random.sample(open_fhs, n_to_close):
        if paired:
            fwd, rev = rand_fh
            fwd.close()
            rev.close()
        else:
            rand_fh.close()


class BarcodeSequenceFastqIterator(collections.abc.Iterable):
    def __init__(self, barcode_generator, sequence_generator,
                 ignore_description_mismatch=False):
        self.barcode_generator = barcode_generator
        self.sequence_generator = sequence_generator
        self.ignore_description_mismatch = ignore_description_mismatch

    def __iter__(self):
        # Adapted from q2-types
        for barcode_record, sequence_record in itertools.zip_longest(
                self.barcode_generator, self.sequence_generator):
            if barcode_record is None:
                raise ValueError('More sequences were provided than barcodes.')
            if sequence_record is None:
                raise ValueError('More barcodes were provided than sequences.')
            # The id or description fields may end with "/read-number", which
            # will differ between the sequence and barcode reads. Confirm that
            # they are identical up until the last /
            barcode_header = _record_to_fastq_header(barcode_record)
            sequence_header = _record_to_fastq_header(sequence_record)

            # confirm that the id fields are equal
            if _trim_id(barcode_header.id) != \
               _trim_id(sequence_header.id):
                raise ValueError(
                    'Mismatched sequence ids: %s and %s' %
                    (_trim_id(barcode_header.id),
                     _trim_id(sequence_header.id)))

            if not self.ignore_description_mismatch:
                # if a description field is present, confirm that they're equal
                if barcode_header.description is None and \
                   sequence_header.description is None:
                    pass
                elif barcode_header.description is None:
                    raise ValueError(
                        'Barcode header lines do not contain description '
                        'fields but sequence header lines do.')
                elif sequence_header.description is None:
                    raise ValueError(
                        'Sequence header lines do not contain description '
                        'fields but barcode header lines do.')
                elif _trim_description(barcode_header.description) != \
                        _trim_description(sequence_header.description):
                    raise ValueError(
                        'Mismatched sequence descriptions: %s and %s' %
                        (_trim_description(barcode_header.description),
                         _trim_description(sequence_header.description)))

            yield barcode_record, sequence_record


class BarcodePairedSequenceFastqIterator(collections.abc.Iterable):
    def __init__(self, barcode_generator, forward_generator,
                 reverse_generator, ignore_description_mismatch=False):
        self.barcode_generator = barcode_generator
        self.forward_generator = forward_generator
        self.reverse_generator = reverse_generator
        self.ignore_description_mismatch = ignore_description_mismatch

    def __iter__(self):
        # Adapted from q2-types
        for barcode_record, forward_record, reverse_record \
                in itertools.zip_longest(self.barcode_generator,
                                         self.forward_generator,
                                         self.reverse_generator):
            if barcode_record is None:
                raise ValueError('More sequences were provided than barcodes.')
            if forward_record is None:
                raise ValueError('More barcodes were provided than '
                                 'forward-sequences.')
            elif reverse_record is None:
                raise ValueError('More barcodes were provided than '
                                 'reverse-sequences.')
            # The id or description fields may end with "/read-number", which
            # will differ between the sequence and barcode reads. Confirm that
            # they are identical up until the last /
            barcode_header = _record_to_fastq_header(barcode_record)
            forward_header = _record_to_fastq_header(forward_record)
            reverse_header = _record_to_fastq_header(reverse_record)

            # confirm that the id fields are equal
            if not (_trim_id(barcode_header.id) ==
                    _trim_id(forward_header.id) ==
                    _trim_id(reverse_header.id)):
                raise ValueError(
                    'Mismatched sequence ids: %s, %s, and %s' %
                    (_trim_id(barcode_header.id),
                     _trim_id(forward_header.id),
                     _trim_id(reverse_header.id)))

            if not self.ignore_description_mismatch:
                # if a description field is present, confirm that they're equal
                if barcode_header.description is None and \
                   forward_header.description is None and \
                   reverse_header.description is None:
                    pass
                elif barcode_header.description is None:
                    raise ValueError(
                        'Barcode header lines do not contain description '
                        'fields but sequence header lines do.')
                elif forward_header.description is None:
                    raise ValueError(
                        'Forward-read header lines do not contain description '
                        'fields but barcode header lines do.')
                elif reverse_header.description is None:
                    raise ValueError(
                        'Reverse-read header lines do not contain description '
                        'fields but barcode header lines do.')
                elif not (_trim_description(barcode_header.description) ==
                          _trim_description(forward_header.description) ==
                          _trim_description(reverse_header.description)):
                    raise ValueError(
                        'Mismatched sequence descriptions: %s, %s, and %s' %
                        (_trim_description(barcode_header.description),
                         _trim_description(forward_header.description),
                         _trim_description(reverse_header.description)))

            yield barcode_record, forward_record, reverse_record


def _make_barcode_map(barcodes, rev_comp_mapping_barcodes):
    barcode_map = {}
    barcode_len = None
    for sample_id, barcode in barcodes.to_series().iteritems():
        if barcode_len is None:
            barcode_len = len(barcode)
        elif len(barcode) != barcode_len:
            raise ValueError('Barcodes of different lengths were detected: '
                             '%d != %d. Variable length barcodes are not '
                             'supported.' % (len(barcode), barcode_len))
        try:
            skbio.DNA(barcode)
        except ValueError as ve:
            if re.match('^ValueError[(]"Invalid characters in sequence[.,'
                        ' \n]*',
                        ve.__repr__()):
                raise ValueError("Invalid characters found in specified "
                                 "barcodes column within metadata file. "
                                 "Please confirm that the column: '%s' "
                                 "contains your per-sample barcodes."
                                 % (barcodes.name))
            else:
                raise

        if rev_comp_mapping_barcodes:
            barcode = str(skbio.DNA(barcode).reverse_complement())

        if barcode in barcode_map:
            raise ValueError('A duplicate barcode was detected. The barcode '
                             '%s was observed for samples %s and %s.'
                             % (barcode, sample_id, barcode_map[barcode]))
        barcode_map[barcode] = sample_id

    return barcode_map, barcode_len


def _write_metadata_yaml(dir_fmt):
    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    dir_fmt.metadata.write_data(metadata, YamlFormat)


def emp_single(seqs: BarcodeSequenceFastqIterator,
               barcodes: qiime2.CategoricalMetadataColumn,
               golay_error_correction: bool = True,
               rev_comp_barcodes: bool = False,
               rev_comp_mapping_barcodes: bool = False,
               ignore_description_mismatch: bool = False
               ) -> (SingleLanePerSampleSingleEndFastqDirFmt,
                     ErrorCorrectionDetailsFmt):
    seqs.ignore_description_mismatch = ignore_description_mismatch
    result = SingleLanePerSampleSingleEndFastqDirFmt()
    barcode_map, barcode_len = _make_barcode_map(
        barcodes, rev_comp_mapping_barcodes)

    if golay_error_correction:
        decoder = GolayDecoder()

    manifest = FastqManifestFormat()
    manifest_fh = manifest.open()
    manifest_fh.write('sample-id,filename,direction\n')
    manifest_fh.write('# direction is not meaningful in this file as these\n')
    manifest_fh.write('# data may be derived from forward, reverse, or \n')
    manifest_fh.write('# joined reads\n')

    per_sample_fastqs = {}

    ec_details_fmt = ErrorCorrectionDetailsFmt()
    ec_details = ECDetails(ec_details_fmt)

    for i, (barcode_record, sequence_record) in enumerate(seqs, start=1):
        raw_barcode_read = barcode_record[1][:barcode_len]
        if rev_comp_barcodes:
            barcode_as_DNA = skbio.DNA(raw_barcode_read)
            raw_barcode_read = str(barcode_as_DNA.reverse_complement())

        if golay_error_correction:
            # A three bit filter is implicitly used by the decoder. See Hamady
            # and Knight 2009 Genome Research for the justification:
            #
            # https://genome.cshlp.org/content/19/7/1141.full
            #
            # Specifically that "...Golay codes of 12 bases can correct all
            # triple-bit errors and detect all quadruple-bit errors."
            barcode_read, ecc_errors = decoder.decode(raw_barcode_read)
            golay_stats = [barcode_read, ecc_errors]
        else:
            barcode_read = raw_barcode_read
            golay_stats = [None, None]

        sample_id = barcode_map.get(barcode_read)

        record = [
            f'record-{i}',
            sample_id,
            barcode_record[0],
            raw_barcode_read,
        ]
        ec_details.write(record + golay_stats)

        if sample_id is None:
            continue

        if sample_id not in per_sample_fastqs:
            # The barcode id, lane number and read number are not relevant
            # here. We might ultimately want to use a dir format other than
            # SingleLanePerSampleSingleEndFastqDirFmt which doesn't care
            # about this information. Similarly, the direction of the read
            # isn't relevant here anymore.
            barcode_id = len(per_sample_fastqs) + 1
            path = result.sequences.path_maker(sample_id=sample_id,
                                               barcode_id=barcode_id,
                                               lane_number=1,
                                               read_number=1)
            _maintain_open_fh_count(per_sample_fastqs)
            per_sample_fastqs[sample_id] = gzip.open(str(path), mode='a')
            manifest_fh.write('%s,%s,%s\n' % (sample_id, path.name, 'forward'))

        if per_sample_fastqs[sample_id].closed:
            _maintain_open_fh_count(per_sample_fastqs)
            per_sample_fastqs[sample_id] = gzip.open(
                per_sample_fastqs[sample_id].name, mode='a')

        fastq_lines = '\n'.join(sequence_record) + '\n'
        fastq_lines = fastq_lines.encode('utf-8')
        per_sample_fastqs[sample_id].write(fastq_lines)

    if len(per_sample_fastqs) == 0:
        raise ValueError('No sequences were mapped to samples. Check that '
                         'your barcodes are in the correct orientation (see '
                         'the rev_comp_barcodes and/or '
                         'rev_comp_mapping_barcodes options). If barcodes are '
                         'NOT Golay format set golay_error_correction '
                         'to False.')

    for fh in per_sample_fastqs.values():
        fh.close()

    manifest_fh.close()
    result.manifest.write_data(manifest, FastqManifestFormat)

    _write_metadata_yaml(result)

    return result, ec_details_fmt


def emp_paired(seqs: BarcodePairedSequenceFastqIterator,
               barcodes: qiime2.CategoricalMetadataColumn,
               golay_error_correction: bool = True,
               rev_comp_barcodes: bool = False,
               rev_comp_mapping_barcodes: bool = False,
               ignore_description_mismatch: bool = False
               ) -> (SingleLanePerSamplePairedEndFastqDirFmt,
                     ErrorCorrectionDetailsFmt):
    seqs.ignore_description_mismatch = ignore_description_mismatch
    result = SingleLanePerSamplePairedEndFastqDirFmt()
    barcode_map, barcode_len = _make_barcode_map(
        barcodes, rev_comp_mapping_barcodes)

    if golay_error_correction:
        decoder = GolayDecoder()

    manifest = FastqManifestFormat()
    manifest_fh = manifest.open()
    manifest_fh.write('sample-id,filename,direction\n')

    per_sample_fastqs = {}

    ec_details_fmt = ErrorCorrectionDetailsFmt()
    ec_details = ECDetails(ec_details_fmt)

    for i, record in enumerate(seqs, start=1):
        barcode_record, forward_record, reverse_record = record
        raw_barcode_read = barcode_record[1][:barcode_len]
        if rev_comp_barcodes:
            barcode_as_DNA = skbio.DNA(raw_barcode_read)
            raw_barcode_read = str(barcode_as_DNA.reverse_complement())

        if golay_error_correction:
            # A three bit filter is implicitly used by the decoder. See Hamady
            # and Knight 2009 Genome Research for the justification:
            #
            # https://genome.cshlp.org/content/19/7/1141.full
            #
            # Specifically that "...Golay codes of 12 bases can correct all
            # triple-bit errors and detect all quadruple-bit errors."
            barcode_read, ecc_errors = decoder.decode(raw_barcode_read)
            golay_stats = [barcode_read, ecc_errors]
        else:
            barcode_read = raw_barcode_read
            golay_stats = [None, None]

        sample_id = barcode_map.get(barcode_read)

        record = [
            f'record-{i}',
            sample_id,
            barcode_record[0],
            raw_barcode_read,
        ]
        ec_details.write(record + golay_stats)

        if sample_id is None:
            continue

        if sample_id not in per_sample_fastqs:
            barcode_id = len(per_sample_fastqs) + 1
            fwd_path = result.sequences.path_maker(sample_id=sample_id,
                                                   barcode_id=barcode_id,
                                                   lane_number=1,
                                                   read_number=1)
            rev_path = result.sequences.path_maker(sample_id=sample_id,
                                                   barcode_id=barcode_id,
                                                   lane_number=1,
                                                   read_number=2)

            _maintain_open_fh_count(per_sample_fastqs, paired=True)
            per_sample_fastqs[sample_id] = (
                gzip.open(str(fwd_path), mode='a'),
                gzip.open(str(rev_path), mode='a')
            )
            manifest_fh.write('%s,%s,%s\n' % (sample_id, fwd_path.name,
                                              'forward'))
            manifest_fh.write('%s,%s,%s\n' % (sample_id, rev_path.name,
                                              'reverse'))

        if per_sample_fastqs[sample_id][0].closed:
            _maintain_open_fh_count(per_sample_fastqs, paired=True)
            fwd, rev = per_sample_fastqs[sample_id]
            per_sample_fastqs[sample_id] = (
                gzip.open(fwd.name, mode='a'),
                gzip.open(rev.name, mode='a')
            )

        fwd, rev = per_sample_fastqs[sample_id]
        fwd.write(('\n'.join(forward_record) + '\n').encode('utf-8'))
        rev.write(('\n'.join(reverse_record) + '\n').encode('utf-8'))

    if len(per_sample_fastqs) == 0:
        raise ValueError('No sequences were mapped to samples. Check that '
                         'your barcodes are in the correct orientation (see '
                         'the rev_comp_barcodes and/or '
                         'rev_comp_mapping_barcodes options). If barcodes are '
                         'NOT Golay format set golay_error_correction '
                         'to False.')

    for fwd, rev in per_sample_fastqs.values():
        fwd.close()
        rev.close()

    manifest_fh.close()
    result.manifest.write_data(manifest, FastqManifestFormat)

    _write_metadata_yaml(result)

    return result, ec_details_fmt
