# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages


setup(
    name="q2-demux",
    version="2017.2.0.dev0",
    packages=find_packages(),
    install_requires=['qiime2 == 2017.2.*', 'q2-types == 2017.2.*',
                      'q2templates == 2017.2.*', 'numpy', 'pandas',
                      'scikit-bio', 'seaborn', 'psutil'],
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
        'q2_demux.test': ['data/**/*'],
        'q2_demux': ['assets/index.html']
    }
)
