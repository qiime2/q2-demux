# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import (
    Plugin, MetadataColumn, Categorical, Bool, Int, Float, Range
)
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality,
    JoinedSequencesWithQuality)

import q2_demux
from ._type import RawSequences, EMPSingleEndSequences, EMPPairedEndSequences
from ._format import (EMPMultiplexedDirFmt,
                      EMPSingleEndDirFmt, EMPSingleEndCasavaDirFmt,
                      EMPPairedEndDirFmt, EMPPairedEndCasavaDirFmt)


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

plugin.register_semantic_types(
    RawSequences, EMPSingleEndSequences, EMPPairedEndSequences)

plugin.register_formats(EMPMultiplexedDirFmt,
                        EMPSingleEndDirFmt, EMPSingleEndCasavaDirFmt,
                        EMPPairedEndDirFmt, EMPPairedEndCasavaDirFmt)

# TODO: remove when aliasing exists
plugin.register_semantic_type_to_format(
    RawSequences,
    artifact_format=EMPSingleEndDirFmt
)

plugin.register_semantic_type_to_format(
    EMPSingleEndSequences,
    artifact_format=EMPSingleEndDirFmt
)


plugin.register_semantic_type_to_format(
    EMPPairedEndSequences,
    artifact_format=EMPPairedEndDirFmt
)


plugin.methods.register_function(
    function=q2_demux.emp_single,
    # TODO: remove RawSequences by creating an alias to EMPSequences
    inputs={'seqs': (RawSequences |
                     EMPSingleEndSequences |
                     EMPPairedEndSequences)},
    parameters={'barcodes': MetadataColumn[Categorical],
                'rev_comp_barcodes': Bool,
                'rev_comp_mapping_barcodes': Bool},
    outputs=[('per_sample_sequences', SampleData[SequencesWithQuality])],
    input_descriptions={
        'seqs': 'The single-end sequences to be demultiplexed.'
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata column containing the per-sample '
                    'barcodes.',
        'rev_comp_barcodes': 'If provided, the barcode sequence reads will be '
                             'reverse complemented prior to demultiplexing.',
        'rev_comp_mapping_barcodes': 'If provided, the barcode sequences in '
                                     'the sample metadata will be reverse '
                                     'complemented prior to demultiplexing.'
    },
    output_descriptions={
        'per_sample_sequences': 'The resulting demultiplexed sequences.'
    },
    name='Demultiplex sequence data generated with the EMP protocol.',
    description=('Demultiplex sequence data (i.e., map barcode reads to '
                 'sample ids) for data generated with the Earth Microbiome '
                 'Project (EMP) amplicon sequencing protocol. Details about '
                 'this protocol can be found at '
                 'http://www.earthmicrobiome.org/protocols-and-standards/')
)

plugin.methods.register_function(
    function=q2_demux.emp_paired,
    inputs={'seqs': EMPPairedEndSequences},
    parameters={'barcodes': MetadataColumn[Categorical],
                'rev_comp_barcodes': Bool,
                'rev_comp_mapping_barcodes': Bool},
    outputs=[
        ('per_sample_sequences', SampleData[PairedEndSequencesWithQuality])
    ],
    input_descriptions={
        'seqs': 'The paired-end sequences to be demultiplexed.'
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata column containing the per-sample '
                    'barcodes.',
        'rev_comp_barcodes': 'If provided, the barcode sequence reads will be '
                             'reverse complemented prior to demultiplexing.',
        'rev_comp_mapping_barcodes': 'If provided, the barcode sequences in '
                                     'the sample metadata will be reverse '
                                     'complemented prior to demultiplexing.'
    },
    output_descriptions={
        'per_sample_sequences': 'The resulting demultiplexed sequences.'
    },
    name=('Demultiplex paired-end sequence data generated with the EMP '
          'protocol.'),
    description=('Demultiplex paired-end sequence data (i.e., map barcode '
                 'reads to sample ids) for data generated with the Earth '
                 'Microbiome Project (EMP) amplicon sequencing protocol. '
                 'Details about this protocol can be found at '
                 'http://www.earthmicrobiome.org/protocols-and-standards/')
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
                 'selected sequences.')
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

importlib.import_module('q2_demux._transformer')
