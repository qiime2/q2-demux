# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from ._demux import emp_single, emp_paired, summarize


__version__ = pkg_resources.get_distribution('q2-demux').version

__all__ = ['emp_single', 'emp_paired', 'summarize']
