# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
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
                           'data/single_sample_multiple_files/*',
                           'data/subsample_single_end/*',
                           'data/subsample_paired_end/*',
                           ],
        'q2_demux': ['_summarize/assets/*.html',
                     '_summarize/assets/dist/*']
    },
    zip_safe=False,
)
