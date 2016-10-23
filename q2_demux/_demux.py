import os.path
import gzip
import yaml
import itertools
import collections
import pkg_resources

import skbio
import pandas as pd
import seaborn as sns

import qiime
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt, FastqManifestFormat, YamlFormat,
    FastqGzFormat)
import q2templates


TEMPLATES = pkg_resources.resource_filename('q2_demux', 'assets')
FastqHeader = collections.namedtuple('FastqHeader', ['id', 'description'])


def _read_fastq_seqs(filepath):
    # This function is adapted from @jairideout's SO post:
    # http://stackoverflow.com/a/39302117/3424666
    fh = gzip.open(filepath, 'rt')
    for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
        yield (seq_header.strip(), seq.strip(), qual_header.strip(),
               qual.strip())


def _trim_trailing_slash(header):
    return header.rsplit('/', 1)[0]


def _record_to_fastq_header(record):
    tokens = record[0][1:].split(' ', maxsplit=1)
    if len(tokens) == 1:
        id, = tokens
        description = None
    else:
        id, description = tokens

    return FastqHeader(id=id, description=description)


class BarcodeSequenceFastqIterator(collections.Iterator):
    def __init__(self, barcode_generator, sequence_generator):
        self.barcode_generator = barcode_generator
        self.sequence_generator = sequence_generator

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
            if _trim_trailing_slash(barcode_header.id) != \
               _trim_trailing_slash(sequence_header.id):
                raise ValueError(
                    'Mismatched sequence ids: %s and %s' %
                    (_trim_trailing_slash(barcode_header.id),
                     _trim_trailing_slash(sequence_header.id)))

            # if a description field is present, confirm that they're equal
            if barcode_header.description is None and \
               sequence_header.description is None:
                pass
            elif barcode_header.description is None:
                raise ValueError(
                    'Barcode header lines do not contain description fields '
                    'but sequence header lines do.')
            elif sequence_header.description is None:
                raise ValueError(
                    'Sequence header lines do not contain description fields '
                    'but barcode header lines do.')
            elif _trim_trailing_slash(barcode_header.description) != \
                    _trim_trailing_slash(sequence_header.description):
                raise ValueError(
                    'Mismatched sequence descriptions: %s and %s' %
                    (_trim_trailing_slash(barcode_header.description),
                     _trim_trailing_slash(sequence_header.description)))

            yield barcode_record, sequence_record

    def __next__(self):
        return next(self.barcode_generator), next(self.sequence_generator)


def summarize(output_dir: str, data: SingleLanePerSampleSingleEndFastqDirFmt) \
        -> None:
    per_sample_fastqs = list(data.sequences.iter_views(FastqGzFormat))
    per_sample_fastq_counts = {}
    for per_sample_fastq in per_sample_fastqs:
        seqs = _read_fastq_seqs(str(per_sample_fastq[1]))
        count = 0
        for seq in seqs:
            count += 1
        per_sample_fastq_counts[per_sample_fastq[0]] = count

    result = pd.Series(per_sample_fastq_counts)
    result.name = 'Sequence count'
    result.index.name = 'Sample filename'
    result.sort_values(inplace=True, ascending=False)
    result.to_csv(os.path.join(output_dir, 'per-sample-fastq-counts.csv'),
                  header=True, index=True)
    ax = sns.distplot(result, kde=False)
    ax.set_xlabel('Number of sequences')
    ax.set_ylabel('Frequency')
    fig = ax.get_figure()
    fig.savefig(os.path.join(output_dir, 'demultiplex-summary.png'))
    fig.savefig(os.path.join(output_dir, 'demultiplex-summary.pdf'))

    html = result.to_frame().to_html(classes='table table-striped table-hover')
    html = html.replace('border="1"', 'border="0"')
    index = os.path.join(TEMPLATES, 'index.html')
    context = {
        'result_data': {
            'min': result.min(),
            'median': result.median(),
            'mean': result.mean(),
            'max': result.max(),
            'sum': result.sum()
        },
        'result': html
    }
    q2templates.render(index, output_dir, context=context)


def emp(seqs: BarcodeSequenceFastqIterator,
        barcodes: qiime.MetadataCategory,
        rev_comp_barcodes: bool=False,
        rev_comp_mapping_barcodes: bool=False) \
        -> SingleLanePerSampleSingleEndFastqDirFmt:
    barcode_map = {}
    barcode_len = None
    for sample_id, barcode in barcodes.to_series().iteritems():
        if barcode_len is None:
            barcode_len = len(barcode)
        elif len(barcode) != barcode_len:
            raise ValueError('Barcodes of different lengths were detected: '
                             '%d != %d. Variable length barcodes are not '
                             'supported.' % (len(barcode), barcode_len))
        if rev_comp_mapping_barcodes:
            barcode = str(skbio.DNA(barcode).reverse_complement())
        if barcode in barcode_map:
            raise ValueError('A duplicate barcode was detected. The barcode '
                             '%s was observed for samples %s and %s.'
                             % (barcode, sample_id, barcode_map[barcode]))
        barcode_map[barcode] = sample_id

    result = SingleLanePerSampleSingleEndFastqDirFmt()

    manifest = FastqManifestFormat()
    manifest_fh = manifest.open()
    manifest_fh.write('sample-id,filename,direction\n')
    manifest_fh.write('# direction is not meaningful in this file as these\n')
    manifest_fh.write('# data may be derived from forward, reverse, or \n')
    manifest_fh.write('# joined reads\n')

    per_sample_fastqs = {}

    for barcode_record, sequence_record in seqs:
        barcode_read = barcode_record[1]
        if rev_comp_barcodes:
            barcode_read = str(skbio.DNA(barcode_read).reverse_complement())
        barcode_read = barcode_read[:barcode_len]

        try:
            sample_id = barcode_map[barcode_read]
        except KeyError:
            # TODO: this should ultimately be logged, but we don't currently
            # have support for that.
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
            per_sample_fastqs[sample_id] = gzip.open(str(path), mode='w')
            manifest_fh.write('%s,%s,%s\n' % (sample_id, path.name, 'forward'))

        fastq_lines = '\n'.join(sequence_record) + '\n'
        fastq_lines = fastq_lines.encode('utf-8')
        per_sample_fastqs[sample_id].write(fastq_lines)

    if len(per_sample_fastqs) == 0:
        raise ValueError('No sequences were mapped to samples. Check that '
                         'your barcodes are in the correct orientation (see '
                         'rev_comp_barcodes and/or rev_comp_mapping_barcodes '
                         'options).')

    for fh in per_sample_fastqs.values():
        fh.close()

    manifest_fh.close()
    result.manifest.write_data(manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    result.metadata.write_data(metadata, YamlFormat)

    return result
