# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import pkg_resources
import shutil
import random
import math

import pandas as pd
import seaborn as sns
import numpy as np

from q2_demux._demux import _read_fastq_seqs
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_demux', '_summarize')


def _decode_qual_to_phred33(qual_str):
    # this function is adapted from scikit-bio
    qual = np.fromstring(qual_str, dtype=np.uint8) - 33
    return qual


# TODO: Remove _PlotQualView once QIIME 2 #220 completed
class _PlotQualView:
    """
    A very simple pass-through view which is made up of a single-end or
    paired-end directory format with a bool indicating if single or paired.
    """
    def __init__(self, directory_format, paired):
        self.directory_format = directory_format
        self.paired = paired


def _link_sample_n_to_file(files, counts, subsample_ns):
    results = collections.defaultdict(list)
    for num in subsample_ns:
        total = 0
        for file in files:
            sample_name = os.path.basename(file).split('_', 1)[0]
            total += counts[sample_name]
            if num < total:
                idx = counts[sample_name] - (total - num)
                results[file].append(idx)
                break
    return results

def _subsample_paired(fastq_map, delta=1):
    qual_sample = collections.defaultdict(list)
    min_seq_len = math.inf
    max_seq_len = 0
    warning = ''
    for fwd, rev, index in fastq_map:
        file_pair = zip(_read_fastq_seqs(fwd), _read_fastq_seqs(rev))
        for i, (fseq, rseq) in enumerate(file_pair):
            min_seqs = min(len(fseq[1]), len(rseq[1]))
            max_seqs = max(len(fseq[1]), len(rseq[1]))
            min_seq_len = min(min_seq_len, min_seqs)
            max_seq_len = max(max_seq_len, max_seqs)
            if i == index[0]:
                qual_sample['forward'].append(_decode_qual_to_phred33(fseq[3]))
                qual_sample['reverse'].append(_decode_qual_to_phred33(rseq[3]))
                index.pop(0)
                if len(index) == 0:
                    break
    if abs(max_seq_len - min_seq_len) > delta:
        warning = (inconsistent_length_template %
            (min_seq_len, max_seq_len))
        for key in qual_sample.keys():
            array_list = qual_sample[key]
            for i in range(len(array_list)):
                array_list[i] = array_list[i][0:min_seq_len:1]
    return [qual_sample, warning]


def _subsample_single(fastq_map, delta=1):
    qual_sample = collections.defaultdict(list)
    min_seq_len = math.inf
    max_seq_len = 0
    warning = ''
    for file, index in fastq_map:
        for i, seq in enumerate(_read_fastq_seqs(file)):
            min_seq_len = min(min_seq_len, len(seq[1]))
            max_seq_len = max(max_seq_len, len(seq[1]))
            if i == index[0]:
                qual_sample['forward'].append(_decode_qual_to_phred33(seq[3]))
                index.pop(0)
                if len(index) == 0:
                    break
    if abs(max_seq_len - min_seq_len) > delta:
        warning = (inconsistent_length_template %
            (min_seq_len, max_seq_len))
        for key in qual_sample.keys():
            array_list = qual_sample[key]
            for i in range(len(array_list)):
                array_list[i] = array_list[i][0:min_seq_len:1]
    return [qual_sample, warning]


def _compute_stats_of_df(df):
    df_stats = df.describe(
        percentiles=[0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98])
    drop_cols = df_stats.index.isin(['count', 'std', 'mean', 'min', 'max'])
    df_stats = df_stats[~drop_cols]
    return df_stats


def summarize(output_dir: str, data: _PlotQualView, n: int=10000, delta=1) -> None:
    paired = data.paired
    data = data.directory_format
    dangers = []
    warnings = []

    manifest = pd.read_csv(os.path.join(str(data), data.manifest.pathspec),
                           header=0, comment='#')
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(data), x))

    fwd = manifest[manifest.direction == 'forward'].filename.tolist()
    rev = manifest[manifest.direction == 'reverse'].filename.tolist()

    per_sample_fastq_counts = {}
    reads = rev if not fwd and rev else fwd
    for file in reads:
        count = 0
        for seq in _read_fastq_seqs(file):
            count += 1
        sample_name = os.path.basename(file).split('_', 1)[0]
        per_sample_fastq_counts[sample_name] = count

    result = pd.Series(per_sample_fastq_counts)
    result.name = 'Sequence count'
    result.index.name = 'Sample name'
    result.sort_values(inplace=True, ascending=False)
    result.to_csv(os.path.join(output_dir, 'per-sample-fastq-counts.csv'),
                  header=True, index=True)
    sequence_count = result.sum()

    if n > sequence_count:
        n = sequence_count
        warnings.append('A subsample value was provided that is greater than '
                        'the amount of sequences across all samples. The plot '
                        'was generated using all available sequences.')

    subsample_ns = sorted(random.sample(range(sequence_count), n))
    link = _link_sample_n_to_file(reads, per_sample_fastq_counts, subsample_ns)
    if paired:
        sample_map = [(file, rev[fwd.index(file)], link[file])
                      for file in link]
        quality_scores = _subsample_paired(sample_map, delta)
    else:
        sample_map = [(file, link[file]) for file in link]
        quality_scores = _subsample_single(sample_map, delta)
    if (quality_scores[1]):
        warnings.append(quality_scores[1])
    quality_scores = quality_scores[0]

    forward_scores = pd.DataFrame(quality_scores['forward'])
    forward_stats = _compute_stats_of_df(forward_scores)

    if (forward_stats.loc['50%'] > 45).any():
        dangers.append('Some of the PHRED quality values are out of range. '
                       'This is likely because an incorrect PHRED offset '
                       'was chosen on import of your raw data. You can learn '
                       'how to choose your PHRED offset during import in the '
                       'importing tutorial.')
    if paired:
        reverse_scores = pd.DataFrame(quality_scores['reverse'])
        reverse_stats = _compute_stats_of_df(reverse_scores)

    show_plot = len(fwd) > 1
    if show_plot:
        ax = sns.distplot(result, kde=False)
        ax.set_xlabel('Number of sequences')
        ax.set_ylabel('Frequency')
        fig = ax.get_figure()
        fig.savefig(os.path.join(output_dir, 'demultiplex-summary.png'))
        fig.savefig(os.path.join(output_dir, 'demultiplex-summary.pdf'))

    html = result.to_frame().to_html(classes='table table-striped table-hover')
    html = html.replace('border="1"', 'border="0"')
    index = os.path.join(TEMPLATES, 'assets', 'index.html')
    overview_template = os.path.join(TEMPLATES, 'assets', 'overview.html')
    quality_template = os.path.join(TEMPLATES, 'assets', 'quality-plot.html')
    context = {
        'result_data': {
            'min': result.min(),
            'median': result.median(),
            'mean': result.mean(),
            'max': result.max(),
            'sum': sequence_count
        },
        'result': html,
        'show_plot': show_plot,
        'paired': paired,
        'tabs': [{'title': 'Overview',
                  'url': 'overview.html'},
                 {'title': 'Interactive Quality Plot',
                  'url': 'quality-plot.html'}],
        'dangers': dangers,
        'warnings': warnings,
        'sample_n': n
    }
    templates = [index, overview_template, quality_template]
    q2templates.render(templates, output_dir, context=context)

    shutil.copytree(os.path.join(TEMPLATES, 'assets', 'app'),
                    os.path.join(output_dir, 'app'))

    with open(os.path.join(output_dir, 'data.jsonp'), 'w') as fh:
        fh.write("app.init(")
        forward_stats.to_json(fh)
        if paired:
            fh.write(',')
            reverse_stats.to_json(fh)
        fh.write(');')


inconsistent_length_template = ('Observed sequences of length %d '
                                'and %d while generating a random '
                                'subsample of sequences. Inconsistent '
                                'length sequences may decrease the '
                                'accuracy of the visualization.')
