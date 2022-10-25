# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


emp_seq_url = 'https://docs.qiime2.org/{epoch}/data/tutorials/' \
              'moving-pictures/emp-single-end-sequences.qza'

demux_url = 'https://docs.qiime2.org/{epoch}/data/tutorials/' \
            'moving-pictures/demux.qza'

metadata_url = 'https://data.qiime2.org/{epoch}/tutorials/' \
               'moving-pictures/sample_metadata.tsv'


def emp_single(use):
    sequences = use.init_artifact_from_url('sequences', emp_seq_url)
    metadata = use.init_metadata_from_url('sample_metadata', metadata_url,
                                          'barcode-sequence')

    demux, correction_details = use.action(
        use.UsageAction('demux', 'emp_single'),
        use.UsageInputs(
            seqs=sequences,
            barcodes=metadata
            ),
        use.UsageOutputNames(
            per_sample_sequences='demux',
            error_correction_details='demux-details'
        )
    )

    demux.assert_output_type('SampleData[SequencesWithQuality]')
    correction_details.assert_output_type('ErrorCorrectionDetails')


def summarize(use):
    demux = use.init_artifact_from_url('demux', demux_url)

    viz, = use.action(
        use.UsageAction('demux', 'summarize'),
        use.UsageInputs(
            data=demux
            ),
        use.UsageOutputNames(
            visualization='visualization'
        )
    )

    viz.assert_output_type('Visualization')
