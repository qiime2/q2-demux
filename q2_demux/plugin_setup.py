# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

import qiime2
import qiime2.plugin
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality)

import q2_demux
from ._type import RawSequences, EMPSingleEndSequences, EMPPairedEndSequences
from ._format import (EMPMultiplexedDirFmt,
                      EMPSingleEndDirFmt, EMPSingleEndCasavaDirFmt,
                      EMPPairedEndDirFmt, EMPPairedEndCasavaDirFmt)


plugin = qiime2.plugin.Plugin(
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
    inputs={'seqs': RawSequences | EMPSingleEndSequences},
    parameters={'barcodes': qiime2.plugin.MetadataCategory,
                'rev_comp_barcodes': qiime2.plugin.Bool,
                'rev_comp_mapping_barcodes': qiime2.plugin.Bool},
    outputs=[('per_sample_sequences', SampleData[SequencesWithQuality])],
    input_descriptions={
        'seqs': 'The single-end sequences to be demultiplexed.'
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata category listing the per-sample '
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
    parameters={'barcodes': qiime2.plugin.MetadataCategory,
                'rev_comp_barcodes': qiime2.plugin.Bool,
                'rev_comp_mapping_barcodes': qiime2.plugin.Bool},
    outputs=[
        ('per_sample_sequences', SampleData[PairedEndSequencesWithQuality])
    ],
    input_descriptions={
        'seqs': 'The paired-end sequences to be demultiplexed.'
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata category listing the per-sample '
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
            SampleData[SequencesWithQuality | PairedEndSequencesWithQuality]},
    parameters={
      'n': qiime2.plugin.Int,
      'delta': qiime2.plugin.Int
    },
    input_descriptions={
        'data': 'The demultiplexed sequences to be summarized.'
    },
    parameter_descriptions={
        'n': ('The number of sequences that should be selected at random for '
              'quality score plots. The quality plots will present the '
              'average positional qualities across all of the sequences '
              'selected. If input sequences are paired end, plots will be '
              'generated for both forward and reverse reads for the same `n` '
              'sequences.'),
        'delta': ('The acceptable delta cutoff for disparities in sequence '
                  'length')
    },
    name='Summarize counts per sample.',
    description=('Summarize counts per sample for all samples, and generate '
                 'interactive positional quality plots based on `n` randomly '
                 'selected sequences.')
)

importlib.import_module('q2_demux._transformer')
