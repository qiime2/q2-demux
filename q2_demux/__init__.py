# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._demux import emp_single, emp_paired
from ._subsample import subsample_single, subsample_paired
from ._summarize import summarize
from ._filter import filter_samples
from ._version import get_versions


__version__ = get_versions()['version']
del get_versions

__all__ = ['emp_single', 'emp_paired', 'summarize',
           'subsample_single', 'subsample_paired', 'filter_samples']
