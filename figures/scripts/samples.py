import collections

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests

from nucperiod import plot as nuceriod_plt, increase


def _sign_SNR(row):
    if row['cross_validation_max'] < 0:
        snr = -np.log2(row['snr'])
    else:
        snr = np.log(row['snr'])
    return snr


def _load(files):

    df = increase.load_d(files)

    df['increase_in'] = df['observed_in'] - df['expected_in']
    df['prop_increase_in'] = df['increase_in'] / df['expected_in']

    df = df.sort_values(by='prop_increase_in', ascending=True)

    df = df[df['nmuts_whole_nucleosome'] > 1000]

    # add Q-value
    df['empirical_pvalue_snr'] = df['empirical_pvalue_snr'].fillna(1)
    df['empirical_pvalue_snr'].replace(0, 0.001, inplace=True)

    qvals = multipletests(df['empirical_pvalue_snr'].tolist(), method='fdr_bh')
    df['qvals_snr'] = qvals[1]

    df['corrected_snr'] = df.apply(_sign_SNR, axis=1)

    return df


def sigmoid(files, tracksheets, ttypes, colors):
    nuceriod_plt.config_params(12)

    df = _load(files)
    toplot = df.merge(tracksheets, how='left')

    toplot['ttype'] = toplot['tumor_name'].map(ttypes)

    ttype_inc = collections.defaultdict(float)
    total_count = 0
    ttype_vals = collections.defaultdict(lambda: collections.defaultdict(list))
    prop_significant = collections.defaultdict(float)
    for ttype, data in toplot.groupby(by='ttype'):
        if len(data) > 10:
            significant = 0
            d = data.sort_values(by='prop_increase_in', ascending=True)
            mean_d = d['prop_increase_in'].median()
            ttype_inc[ttype] = mean_d
            for i, row in d.iterrows():
                ttype_vals[ttype]['Prop'].append(row['prop_increase_in'])
                total_count += 1
                if (row['qvals_snr'] < 0.1) & (row['snr'] > 8):
                    significant += 1
                    c = 'darkred'
                else:
                    c = 'grey'
                ttype_vals[ttype]['col'].append(c)
            prop_significant[ttype] = 100 * significant / len(data)

    fig, ax = plt.subplots(nrows=1, ncols=len(ttype_inc), figsize=(17, 3.5), sharey=True)

    ax[0].set_ylabel('Relative increase in mutation rate')
    ax[0].yaxis.set_ticks(np.arange(-0.3, 0.3, 0.1))

    for index, t in enumerate(sorted(ttype_inc, key=ttype_inc.get)):
        count = 0

        for ix, val in enumerate(ttype_vals[t]['Prop']):
            alpha = 0.3
            if ttype_vals[t]['col'][ix] == 'darkred':
                alpha = 1
            ax[index].scatter(count, val, color=ttype_vals[t]['col'][ix], s=10, lw=0, alpha=alpha)
            count += 1
        ax[index].hlines(ttype_inc[t], count / 2 - count * 0.4 / 2, count / 2 + count * 0.4 / 2, lw=2,
                         color='darkgreen')
        ax[index].spines['top'].set_visible(False)
        ax[index].spines['right'].set_visible(False)
        ax[index].spines['bottom'].set_visible(False)
        ax[index].spines['left'].set_visible(False)
        ax[index].set_xlabel(t, rotation=90)
        ax[index].xaxis.set_ticks_position('none')

        if index > 0:
            ax[index].yaxis.set_ticks_position('none')

        labels = [item.get_text() for item in ax[index].get_xticklabels()]
        empty_string_labels = [''] * len(labels)
        ax[index].set_xticklabels(empty_string_labels)
        ax[index].text(ttype_inc[t], 0.22, 'n={}\n{}%'.format(count, round(prop_significant[t], 1)))
        ax[index].set_ylim(-0.2, 0.3)

        ax[index].add_patch(plt.Rectangle((0, -0.8), count, 0.03, color=colors[t], lw=1, clip_on=False, linewidth=0))

    red_dot = mlines.Line2D([], [], color='darkred', marker='o', linestyle='None', markersize=5,
                            label='Significant sample')

    grey_dot = mlines.Line2D([], [], color='grey', marker='o', linestyle='None', markersize=5,
                             label='Non significant sample')

    median_line = mlines.Line2D([], [], color='darkgreen', marker='_', linestyle='None', lw=40, markersize=10,
                                label='Median')

    plt.legend(handles=[red_dot, grey_dot, median_line], bbox_to_anchor=[1.1, 1.1])


def scatter(files, tracksheets, ttypes, colors):
    nuceriod_plt.config_params(12)

    df = _load(files)
    toplot = df.merge(tracksheets, how='left')

    toplot['ttype'] = toplot['tumor_name'].map(ttypes)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 4))

    ax.set_xscale('log')

    for i, row in toplot[(toplot['qvals_snr'] < 0.1) & (toplot['snr'] > 8)].iterrows():
        ax.scatter(row['muts'], row['prop_increase_in'], s=14, color=colors[row['ttype']])

    plt.title('Significant')
    plt.ylabel('Proportion of increase minor in')
    plt.xlabel('Number of mutations (log)')
    plt.hlines(0, 0, 1000000, linestyle='--', color='grey', alpha=0.6)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylim(-0.25, 0.25)


def _get_project(name):
    if '-' in name:
        return name.split('-')[1]
    else:
        return '505'


def generate_table(files, tracksheets, ttypes):

    df = _load(files)
    toplot = df.merge(tracksheets, how='left')

    toplot['ttype'] = toplot['tumor_name'].map(ttypes)

    toplot['project'] = toplot['tumor_name'].apply(_get_project)

    table_out = toplot[['name', 'ttype', 'tumor_name', 'project', 'peak', 'snr',
                        'empirical_pvalue_snr', 'qvals_snr', 'cross_validation_max']].copy()

    table_out.columns = ['Sample', 'Tumor Name', 'Cohort', 'Project', 'Peak', 'SNR', 'P-value', 'Q-value', 'Phase']

    return table_out.sort_values(by='Sample')


def _load_contributions(files):
    to_concat = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        to_concat.append(df)
    fin = pd.concat(to_concat, ignore_index=True)
    sigicgc = fin.drop(['SSE', 'mutation_count'], axis=1).set_index('sample_id')
    d = sigicgc.to_dict(orient='index')
    return d


def piechart(files, colors_d, samples):

    signatures_contribution = _load_contributions(files)

    fig, ax = plt.subplots(nrows=1, ncols=len(samples), figsize=(10, 2))

    for ix, s in enumerate(samples):
        prop_samples = signatures_contribution[s]
        labels = list(prop_samples.keys())
        sizes = list(prop_samples.values())
        tot = sum([float(i) for i in sizes])
        sizes_real = [float(s) / tot for s in sizes]
        colors = [colors_d[i.replace('.', '_')] for i in labels]

        ax[ix].pie(sizes_real, colors=colors, shadow=True, startangle=90)
        ax[ix].axis('equal')
