# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType

# TODO: migrate these to q2-types someday
RawSequences = SemanticType('RawSequences')

EMPSingleEndSequences = SemanticType('EMPSingleEndSequences')

EMPPairedEndSequences = SemanticType('EMPPairedEndSequences')
