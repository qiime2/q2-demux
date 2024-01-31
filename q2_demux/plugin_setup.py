# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime2.plugin import (
    Plugin, Metadata, MetadataColumn, Categorical, Bool, Str, Int, Float,
    List, Collection, Range, Citations, TypeMatch
)

from q2_types.sample_data import SampleData
from q2_types.metadata import ImmutableMetadata
from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality,
    JoinedSequencesWithQuality)

from q2_types.multiplexed_sequences import (RawSequences,
                                            EMPSingleEndSequences,
                                            EMPPairedEndSequences,
                                            ErrorCorrectionDetails)

import q2_demux
import q2_demux._examples as ex

citations = Citations.load('citations.bib', package='q2_demux')

plugin = Plugin(
    name='demux',
    version=q2_demux.__version__,
    website='https://github.com/qiime2/q2-demux',
    package='q2_demux',
    description=('This QIIME 2 plugin supports demultiplexing of '
                 'single-end and paired-end sequence reads and '
                 'visualization of sequence quality information.'),
    short_description='Plugin for demultiplexing & viewing sequence quality.'
)


plugin.methods.register_function(
    function=q2_demux.emp_single,
    # TODO: remove RawSequences by creating an alias to EMPSequences
    inputs={'seqs': (RawSequences |
                     EMPSingleEndSequences |
                     EMPPairedEndSequences)},
    parameters={'barcodes': MetadataColumn[Categorical],
                'golay_error_correction': Bool,
                'rev_comp_barcodes': Bool,
                'rev_comp_mapping_barcodes': Bool,
                'ignore_description_mismatch': Bool},
    outputs=[('per_sample_sequences', SampleData[SequencesWithQuality]),
             ('error_correction_details', ErrorCorrectionDetails)],
    input_descriptions={
        'seqs': 'The single-end sequences to be demultiplexed.'
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata column containing the per-sample '
                    'barcodes.',
        'golay_error_correction': 'Perform 12nt Golay error correction on the '
                                  'barcode reads.',
        'rev_comp_barcodes': 'If provided, the barcode sequence reads will be '
                             'reverse complemented prior to demultiplexing.',
        'rev_comp_mapping_barcodes': 'If provided, the barcode sequences in '
                                     'the sample metadata will be reverse '
                                     'complemented prior to demultiplexing.',
        'ignore_description_mismatch': 'If enabled, ignore mismatches in '
                                       'sequence record description fields.'
    },
    output_descriptions={
        'per_sample_sequences': 'The resulting demultiplexed sequences.',
        'error_correction_details': 'Detail about the barcode error '
                                    'corrections.'
    },
    name='Demultiplex sequence data generated with the EMP protocol.',
    description=('Demultiplex sequence data (i.e., map barcode reads to '
                 'sample ids) for data generated with the Earth Microbiome '
                 'Project (EMP) amplicon sequencing protocol. Details about '
                 'this protocol can be found at '
                 'http://www.earthmicrobiome.org/protocols-and-standards/'),
    examples={'demux': ex.emp_single},
    citations=[
        citations['hamady2008'],
        citations['hamady2009']]
)

plugin.methods.register_function(
    function=q2_demux.emp_paired,
    inputs={'seqs': EMPPairedEndSequences},
    parameters={'barcodes': MetadataColumn[Categorical],
                'golay_error_correction': Bool,
                'rev_comp_barcodes': Bool,
                'rev_comp_mapping_barcodes': Bool,
                'ignore_description_mismatch': Bool},
    outputs=[
        ('per_sample_sequences', SampleData[PairedEndSequencesWithQuality]),
        ('error_correction_details', ErrorCorrectionDetails),
    ],
    input_descriptions={
        'seqs': 'The paired-end sequences to be demultiplexed.'
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata column containing the per-sample '
                    'barcodes.',
        'golay_error_correction': 'Perform 12nt Golay error correction on the '
                                  'barcode reads.',
        'rev_comp_barcodes': 'If provided, the barcode sequence reads will be '
                             'reverse complemented prior to demultiplexing.',
        'rev_comp_mapping_barcodes': 'If provided, the barcode sequences in '
                                     'the sample metadata will be reverse '
                                     'complemented prior to demultiplexing.',
        'ignore_description_mismatch': 'If enabled, ignore mismatches in '
                                       'sequence record description fields.'
    },
    output_descriptions={
        'per_sample_sequences': 'The resulting demultiplexed sequences.',
        'error_correction_details': 'Detail about the barcode error '
                                    'corrections.'
    },
    name=('Demultiplex paired-end sequence data generated with the EMP '
          'protocol.'),
    description=('Demultiplex paired-end sequence data (i.e., map barcode '
                 'reads to sample ids) for data generated with the Earth '
                 'Microbiome Project (EMP) amplicon sequencing protocol. '
                 'Details about this protocol can be found at '
                 'http://www.earthmicrobiome.org/protocols-and-standards/'),
    citations=[
        citations['hamady2008'],
        citations['hamady2009']]
)

demux_description = 'The demultiplexed sequences to partition.'
num_partitions_description = 'The number of partitions to split the' \
                             ' demultiplexed sequences into. Defaults to' \
                             ' partitioning into individual samples.'
partitioned_demux_description = 'The partitioned demultiplexed sequences.'

plugin.methods.register_function(
    function=q2_demux.partition_samples_single,
    inputs={'demux': SampleData[SequencesWithQuality]},
    parameters={'num_partitions': Int % Range(1, None)},
    outputs=[
        ('partitioned_demux', Collection[SampleData[SequencesWithQuality]]),
    ],
    input_descriptions={
        'demux': demux_description
    },
    parameter_descriptions={
        'num_partitions': num_partitions_description
    },
    output_descriptions={
        'partitioned_demux': partitioned_demux_description
    },
    name='Split demultiplexed sequence data into partitions.',
    description=('Partition demultiplexed single end sequences into '
                 'individual samples or the number of partitions specified.'),
)

plugin.methods.register_function(
    function=q2_demux.partition_samples_paired,
    inputs={'demux': SampleData[PairedEndSequencesWithQuality]},
    parameters={'num_partitions': Int % Range(1, None)},
    outputs=[
        ('partitioned_demux',
         Collection[SampleData[PairedEndSequencesWithQuality]]),
    ],
    input_descriptions={
        'demux': demux_description
    },
    parameter_descriptions={
        'num_partitions': num_partitions_description
    },
    output_descriptions={
        'partitioned_demux': partitioned_demux_description
    },
    name='Split demultiplexed sequence data into partitions.',
    description=('Partition demultiplexed paired end sequences into '
                 'individual samples or the number of partitions specified.'),
)

plugin.visualizers.register_function(
    function=q2_demux.summarize,
    inputs={'data':
            SampleData[SequencesWithQuality |
                       PairedEndSequencesWithQuality |
                       JoinedSequencesWithQuality]},
    parameters={'n': Int},
    input_descriptions={
        'data': 'The demultiplexed sequences to be summarized.'
    },
    parameter_descriptions={
        'n': ('The number of sequences that should be selected at random for '
              'quality score plots. The quality plots will present the '
              'average positional qualities across all of the sequences '
              'selected. If input sequences are paired end, plots will be '
              'generated for both forward and reverse reads for the same `n` '
              'sequences.')
    },
    name='Summarize counts per sample.',
    description=('Summarize counts per sample for all samples, and generate '
                 'interactive positional quality plots based on `n` randomly '
                 'selected sequences.'),
    examples={'demux': ex.summarize}
)

plugin.methods.register_function(
    function=q2_demux.tabulate_read_counts,
    inputs={'sequences':
            List[SampleData[SequencesWithQuality |
                            PairedEndSequencesWithQuality |
                            JoinedSequencesWithQuality]]},
    parameters={},
    outputs=[
        ('counts', ImmutableMetadata)
    ],
    input_descriptions={
        'sequences': 'One or more collections of demultiplexed sequences.'
    },
    parameter_descriptions={},
    name='Tabulate counts per sample',
    description=('Generate a per-sample count of sequence reads.'),
    examples={}
)

plugin.methods.register_function(
    function=q2_demux.subsample_single,
    inputs={'sequences': SampleData[SequencesWithQuality |
                                    PairedEndSequencesWithQuality]},
    parameters={'fraction': Float % Range(0, 1,
                                          inclusive_start=False,
                                          inclusive_end=False)},
    outputs=[
        ('subsampled_sequences', SampleData[SequencesWithQuality])
    ],
    input_descriptions={
        'sequences': 'The demultiplexed sequences to be subsampled.'
    },
    parameter_descriptions={
        'fraction': ('The fraction of sequences to retain in subsample.')
    },
    output_descriptions={
        'subsampled_sequences': 'The subsampled sequences.'
    },
    name='Subsample single-end sequences without replacement.',
    description=('Generate a random subsample of single-end sequences '
                 'containing approximately the fraction of input sequences '
                 'specified by the fraction parameter. The number of output '
                 'samples will always be equal to the number of input '
                 'samples, even if some of those samples contain no '
                 'sequences after subsampling.')
)

plugin.methods.register_function(
    function=q2_demux.subsample_paired,
    inputs={'sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={'fraction': Float % Range(0, 1,
                                          inclusive_start=False,
                                          inclusive_end=False)},
    outputs=[
        ('subsampled_sequences', SampleData[PairedEndSequencesWithQuality])
    ],
    input_descriptions={
        'sequences': 'The demultiplexed sequences to be subsampled.'
    },
    parameter_descriptions={
        'fraction': ('The fraction of sequences to retain in subsample.')
    },
    output_descriptions={
        'subsampled_sequences': 'The subsampled sequences.'
    },
    name='Subsample paired-end sequences without replacement.',
    description=('Generate a random subsample of paired-end sequences '
                 'containing approximately the fraction of input sequences '
                 'specified by the fraction parameter. The number of output '
                 'samples will always be equal to the number of input '
                 'samples, even if some of those samples contain no '
                 'sequences after subsampling.')
)

T = TypeMatch([SequencesWithQuality, PairedEndSequencesWithQuality,
               JoinedSequencesWithQuality])
plugin.methods.register_function(
    function=q2_demux.filter_samples,
    inputs={'demux': SampleData[T]},
    parameters={'metadata': Metadata,
                'where': Str,
                'exclude_ids': Bool},
    outputs=[
        ('filtered_demux', SampleData[T])
    ],
    input_descriptions={
        'demux': 'The demultiplexed data from which samples should be '
        'filtered.'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata indicating which sample ids to filter. '
                    'The optional `where` parameter may be used to filter ids '
                    'based on specified conditions in the metadata. The '
                    'optional `exclude_ids` parameter may be used to exclude '
                    'the ids specified in the metadata from the filter.',
        'where': 'Optional SQLite WHERE clause specifying sample metadata '
                 'criteria that must be met to be included in the filtered '
                 'data. If not provided, all samples in `metadata` that are '
                 'also in the demultiplexed data will be retained.',
        'exclude_ids': 'Defaults to False. If True, the samples selected by '
                       'the `metadata` and optional `where` parameter will be '
                       'excluded from the filtered data.',
    },
    output_descriptions={
        'filtered_demux': 'Filtered demultiplexed data.'
    },
    name='Filter samples out of demultiplexed data.',
    description='Filter samples indicated in given metadata out of '
                'demultiplexed data. Specific samples can be further selected '
                'with the WHERE clause, and the `exclude_ids` parameter '
                'allows for filtering of all samples not specified.',
)


