# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_demux._demux import (emp_single, emp_paired, partition_samples_single,
                             partition_samples_paired)
from q2_demux._subsample import subsample_single, subsample_paired
from q2_demux._summarize import summarize
from q2_demux._filter import filter_samples
from q2_demux._tabulate import tabulate_read_counts

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = '0.0.0+notfound'

__all__ = ['emp_single', 'emp_paired', 'partition_samples_single',
           'partition_samples_paired', 'summarize', 'subsample_single',
           'subsample_paired', 'filter_samples', 'tabulate_read_counts']
