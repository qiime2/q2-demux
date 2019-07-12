# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer


setup(
    name="q2-demux",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    url="https://qiime2.org",
    license="BSD-3-Clause",
    description="Maps sequence barcodes to samples.",
    entry_points={
        "qiime2.plugins":
        ["q2-demux=q2_demux.plugin_setup:plugin"]
    },
    package_data={
        'q2_demux.tests': ['data/bad/*',
                           'data/emp_multiplexed/*',
                           'data/emp_multiplexed_single_end/*',
                           'data/error_correction_details/*',
                           'data/single_sample_multiple_files/*',
                           'data/subsample_single_end/*',
                           'data/subsample_paired_end/*',
                           'data/filter_samples_single_end/*',
                           'data/filter_samples_single_end/dir_fmt/*',
                           'data/filter_samples_paired_end/*',
                           'data/filter_samples_paired_end/dir_fmt/*',
                           'data/summarize_empty/*'
                           'data/summarize_empty/empty_SampleDataSingleEnd/*',
                           'data/summarize_empty/'
                           'empty_SampleDataPairedEndForward/*',
                           'data/summarize_empty/'
                           'empty_SampleDataPairedEndReverse/*',
                           ],
        'q2_demux': ['_summarize/assets/*.html',
                     '_summarize/assets/dist/*',
                     'citations.bib'],
    },
    zip_safe=False,
)
