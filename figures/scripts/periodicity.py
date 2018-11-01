import json
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats import multitest

from nucperiod import plot as nucperiod_plt, spectral, non_linear


def name_convert(name):
    genus, species = name.split('_', 1)
    species = species.replace('_', ' ')
    return ' '.join([genus.capitalize(), species])


def build_families(folder):
    """Create dict with organisms arranged by families"""

    organisms = defaultdict(dict)
    for k in os.listdir(folder):
        if '_' in k and not k.startswith('.'):  # avoid conflicts with other possible folders as scripts or .ipynb_checkpoints
            family, species = k.split('_', 1)
            species = name_convert(species)
        else:
            continue
        organisms[family][species] = os.path.join(folder, k)
    return organisms


def plot_single(table, title):
    nucperiod_plt.config_params(font_size=10)

    YLIM = 60

    table = pd.read_csv(table, sep='\t')

    fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(20, 12))

    snrs = []
    df = table[table['maxp'] < 19.5]
    df = df[df['maxp'] > 5.5]
    for i in range(6, 20):
        snrs.append(df['snr_{0}p'.format(str(i))].values)
    ax[0, 0].boxplot(snrs)
    ax[0, 0].set_xticklabels(range(6, 20))
    ax[0, 0].set_xlabel('Period (bp)')
    ax[0, 0].set_ylabel('SNR')
    ax[0, 0].set_title('SNR for all 1 Mb chunks')

    snrs = []
    counts = []
    for i in range(6, 20):
        df = table[table['maxp'] < i + 0.5]
        df = df[df['maxp'] > i - 0.5]
        counts.append(len(df))
        snrs.append(df['snr_{0}p'.format(str(i))].values)
    ax[1, 0].boxplot(snrs)
    ax[1, 0].plot(range(1, 15), [v / 10 for v in counts], label='x10 no. chunks')
    ax[1, 0].set_xticklabels(range(6, 20))
    ax[1, 0].set_xlabel('Period (bp)')
    ax[1, 0].set_ylabel('SNR')
    ax[1, 0].set_title('SNR for 1 Mb chunks\n peaking at each period')
    ax[1, 0].legend()

    snrs = []
    df = table[table['maxp'] < 19.5]
    df = df[df['maxp'] > 5.5]
    for i in range(6, 20):
        snrs.append(df['fold_power_increase_{0}p'.format(str(i))].values)
    ax[0, 1].boxplot(snrs)
    ax[0, 1].set_xticklabels(range(6, 20))
    ax[0, 1].set_ylim(-3, YLIM)
    ax[0, 1].set_xlabel('Period (bp)')
    ax[0, 1].set_ylabel('Fold Power Increase')
    ax[0, 1].set_title('Fold power increase for all 1 Mb chunks')

    snrs = []
    for i in range(6, 20):
        df = table[table['maxp'] < i + 0.5]
        df = df[df['maxp'] > i - 0.5]
        snrs.append(df['fold_power_increase_{0}p'.format(str(i))].values)
    ax[1, 1].boxplot(np.array(snrs))
    ax[1, 1].plot(range(1, 15), [v / 10 for v in counts], label='x10 no. chunks')
    ax[1, 1].set_xticklabels(range(6, 20))
    ax[1, 1].set_ylim(-3, YLIM)
    ax[1, 1].set_xlabel('Period (bp)')
    ax[1, 1].set_ylabel('Fold Power Increase')
    ax[1, 1].set_title('Fold Power Increase for 1 Mb chunks\n peaking at each period')
    ax[1, 1].legend()

    snrs = []
    df = table[table['maxp'] < 19.5]
    df = df[df['maxp'] > 5.5]
    for i in range(6, 20):
        snrs.append(df['fold_snr_increase_{0}p'.format(str(i))].values)
    ax[0, 2].boxplot(snrs)
    ax[0, 2].set_xticklabels(range(6, 20))
    ax[0, 2].set_ylim(-3, YLIM)
    ax[0, 2].set_xlabel('Period (bp)')
    ax[0, 2].set_ylabel('Fold SNR Increase')
    ax[0, 2].set_title('Fold SNR Increase for all 1 Mb chunks')

    snrs = []
    for i in range(6, 20):
        df = table[table['maxp'] < i + 0.5]
        df = df[df['maxp'] > i - 0.5]
        snrs.append(df['fold_snr_increase_{0}p'.format(str(i))].values)
    ax[1, 2].boxplot(snrs)
    ax[1, 2].plot(range(1, 15), [v / 10 for v in counts], label='x10 no. chunks')
    ax[1, 2].set_xticklabels(range(6, 20))
    ax[1, 2].set_ylim(-3, YLIM)
    ax[1, 2].set_xlabel('Period (bp)')
    ax[1, 2].set_ylabel('Fold SNR Increase')
    ax[1, 2].set_title('Fold SNR Increase for 1 Mb chunks\n peaking at each period')
    ax[1, 2].legend()

    discoveries = []
    for i in range(6, 20):
        discoveries.append(len(table[table['logpval_power_{0}p'.format(str(i))] == 2]))
    barlist = ax[0, 3].bar(list(range(6, 20)), discoveries)
    barlist[list(range(6, 20)).index(10)].set_color('r')
    ax[0, 3].hlines(np.mean(np.array(discoveries)), 5.6, 19.4, linestyles='dashed', colors='grey')
    ax[0, 3].set_xticks(range(6, 20))
    ax[0, 3].set_xlabel('Period (bp)')
    ax[0, 3].set_ylabel('No. chunks significantly high in power')
    ax[0, 3].set_title('Power Enrichment')

    discoveries = []
    for i in range(6, 20):
        discoveries.append(len(table[table['logpval_snr_{0}p'.format(str(i))] == 2]))
    barlist = ax[1, 3].bar(list(range(6, 20)), discoveries)
    barlist[list(range(6, 20)).index(10)].set_color('r')
    ax[1, 3].hlines(np.mean(np.array(discoveries)), 5.6, 19.4, linestyles='dashed', colors='grey')
    ax[1, 3].set_xticks(range(6, 20))
    ax[1, 3].set_ylabel('No. chunks significantly high in SNR')
    ax[1, 3].set_xlabel('Period (bp)')
    ax[1, 3].set_title('SNR Enrichment')

    fig.suptitle(title)
    plt.rcParams['savefig.facecolor'] = fig.get_facecolor()


def scatter(df, feature1, feature2, selected_organisms=None, **kwargs):
    just_once = True
    if selected_organisms is None:
        selected_organisms = []
    nucperiod_plt.config_params(font_size=24)
    fig, axes = plt.subplots(1, 1, figsize=(12, 14))
    phyllum = ['protists', 'fungi', 'plants', 'vertebrates', 'insects', 'nematodes', 'deuterostomes']
    color_label = dict(zip(phyllum, ['blue', 'grey', 'green', 'pink', 'black', 'cyan', 'orange']))
    for org_type in phyllum:
        ds = df[df['org_type'] == org_type]
        linewidths_normal = [1 for _ in ds[feature1].values]
        linewidths_snr = [2 if a < 1e-2 and b > 0 else 0. for a, b in zip(ds['qval_power_enrichment'].values, ds[feature1].values)]
        if just_once:
            axes.scatter([], [], s=300, linewidths=linewidths_snr,
                         edgecolor='#8b0000ff', color='None', label='q-value < 0.01')
            just_once = False
        axes.scatter(ds[feature1].values, ds[feature2].values, s=350, linewidths=linewidths_snr,
                     edgecolor='#8b0000ff', color='None')
        axes.scatter(ds[feature1].values, ds[feature2].values, s=150, linewidths=linewidths_normal,
                     edgecolor='black', color=color_label[org_type], label=org_type)
    for i, txt in enumerate(df.index.values):
        if txt in selected_organisms:
            axes.annotate(txt, (df.loc[txt, feature1] + 0.02, df.loc[txt, feature2] + 0.02))
    axes.set_ylabel('Proportion of 1 Mb chunks with with MP at {0} $\pm$ 0.5 bp (ratio)'.format(kwargs['period']))
    axes.set_xlabel('Power Enrichment at {0} $\pm$ 0.5 bp (odds ratio)'.format(kwargs['period']))
    # Power enrichment: power significance at {0} $\pm$ 0.5 bp compared to other periods (odds ratio)
    axes.vlines(0, 0, 1, colors='red', linestyles='dashed', alpha=0.5)
    axes.legend(loc=(1.03, 0.712))
    axes.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    axes.set_xticklabels([0.001, 0.01, 0.1, 1, 10, 100, 1000])
    axes.set_title('Period = {0} bp'.format(kwargs['period']))
    if kwargs:
        xmin = kwargs['xmin']
        ymin = kwargs['ymin']
        xmax = kwargs['xmax']
        ymax = kwargs['ymax']
        axes.set_xlim(xmin, xmax)
        axes.set_ylim(ymin, ymax)

    plt.rcParams['savefig.facecolor'] = fig.get_facecolor()


def period_values(families, ref_period):
    df_list = []

    for family in families:
        for species, w in families[family].items():
            table = os.path.join(w, 'metrics.tsv')
            if os.path.exists(table):
                df = pd.read_csv(table, sep='\t')
                df['org_type'] = family
                df['species'] = species
                df = df[df['period'] == ref_period]
                df_list.append(df)

    dh = pd.concat(df_list, axis=0)

    dh['qval_snr'] = dh['pval_snr_all_chunks'].apply(lambda x: multitest.multipletests(x, method='holm')[1][0])
    dh['qval_fold_power_all_chunks'] = dh['pval_fold_power_all_chunks'].apply(
        lambda x: multitest.multipletests(x, method='holm')[1][0])
    dh['qval_fold_snr_all_chunks'] = dh['pval_fold_snr_all_chunks'].apply(
        lambda x: multitest.multipletests(x, method='holm')[1][0])
    dh['qval_fold_power_peaks'] = dh['pval_fold_power_peaks'].apply(
        lambda x: multitest.multipletests(x, method='holm')[1][0])
    dh['qval_fold_snr_peaks'] = dh['pval_fold_snr_peaks'].apply(
        lambda x: multitest.multipletests(x, method='holm')[1][0])
    dh['qval_power_enrichment'] = dh['pval_power_enrichment'].apply(
        lambda x: multitest.multipletests(x, method='holm')[1][0])
    dh['qval_snr_enrichment'] = dh['pval_snr_enrichment'].apply(
        lambda x: multitest.multipletests(x, method='holm')[1][0])
    dh = dh.replace([np.inf, -np.inf], 100)

    return dh


def dtft_spectrum_plot(signal, ax, norm=True, low_t=30, high_t=100, low_p=5, high_p=20, title='Title', alpha=1, color='green', label=None):
    """
    Args:
        signal: numpy array with discrete-time signal values
        ax: matplotlib axes
        low_t: int: lower bound for t in DTFT
        high_t: int: upper bound for t in DTFT
        low_p: float: lower bound for p in DTFT
        high_p: float: upper bound for p in DTFT
        title: str: plot title
    Returns:
        plot in axes
    """
    x, y, snr, peak = spectral.compute(signal, norm=norm, center=10, low_t=low_t, high_t=high_t, low_p=low_p, high_p=high_p)
    ax.plot(x, y, alpha=alpha, linewidth=3., color=color, label=label)
    ax.set_title(title)
    ax.set_xticks(np.arange(low_p, high_p))
    return x, y


def autocorrelation(original_file, simulated_files):
    nucperiod_plt.config_params(font_size=14)

    with open(original_file, 'rt') as f:
        mydict = json.load(f)

    pair_count_array = np.array(mydict['pair_count'])
    motif_count = mydict['motif_count']
    len_chunk = mydict['chunk_len']

    # define the signal
    signal = np.array([(v / len_chunk) / (motif_count / len_chunk) ** 2 for v in pair_count_array])

    # define figure
    figsize = (10, 5)
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.set_title('motif autocorrelation')

    # raw signal
    ax.plot(range(4, len(signal)), signal[4:], alpha=0.5, label='raw')

    # 3-smoothing
    signal = spectral.mean_smooth(signal, 3)
    ax.plot(range(4, len(signal)), signal[4:], linewidth=3.0, label='3-smoothed')

    counter = 0

    # detrended smooth autocorrelation
    initial_values_dict = {'a_0': np.mean(signal[4:]), 'a_1': 0., 'a_2': 0.}
    params, obj_func = non_linear.create_quadratic_model(initial_values_dict)
    x = np.arange(len(signal[4:]))
    non_linear_fitter = non_linear.NonLinearFitting(obj_func, params, x, signal[4:])
    _, predicted = non_linear_fitter.least_squares()

    # with quadratic least-squares fit
    ax.plot(range(4, len(signal)), predicted, '--', label='quadratic trend')

    for file in simulated_files:
        counter += 1

        with open(file, 'rt') as f:
            random_chunk = json.load(f)

        pc = np.array(random_chunk['pair_count'])
        mc = random_chunk['motif_count']
        len_random_chunk = random_chunk['chunk_len']

        random_signal = np.array([(v / len_chunk) / (mc / len_random_chunk) ** 2 for v in pc])
        random_signal = spectral.mean_smooth(random_signal, 3)
        if counter == 1:
            label = 'randomized'
        else:
            label = None
        ax.plot(range(4, len(random_signal)), random_signal[4:], linewidth=3.0,
                alpha=0.3, color='grey', label=label)

    ax.legend()
    plt.rcParams['savefig.facecolor'] = fig.get_facecolor()


def spectrum(original_file, simulated_files):
    nucperiod_plt.config_params(font_size=14)

    with open(original_file, 'rt') as f:
        mydict = json.load(f)

    pair_count_array = np.array(mydict['pair_count'])
    motif_count = mydict['motif_count']
    len_chunk = mydict['chunk_len']

    # define the signal
    signal = np.array([(v / len_chunk) / (motif_count / len_chunk) ** 2 for v in pair_count_array])

    # define figure
    figsize = (10, 5)
    fig_spec, ax_spec = plt.subplots(1, 1, figsize=figsize)
    ax_spec.set_xlabel('period (bp)')
    ax_spec.set_ylabel('power')

    # 3-smoothing
    signal = spectral.mean_smooth(signal, 3)

    counter = 0

    # detrended smooth autocorrelation
    initial_values_dict = {'a_0': np.mean(signal[4:]), 'a_1': 0., 'a_2': 0.}
    params, obj_func = non_linear.create_quadratic_model(initial_values_dict)
    x = np.arange(len(signal[4:]))
    non_linear_fitter = non_linear.NonLinearFitting(obj_func, params, x, signal[4:])
    _, predicted = non_linear_fitter.least_squares()
    signal_detrended = signal[4:] - predicted

    for file in simulated_files:
        counter += 1

        with open(file, 'rt') as f:
            random_chunk = json.load(f)

        pc = np.array(random_chunk['pair_count'])
        mc = random_chunk['motif_count']
        len_random_chunk = random_chunk['chunk_len']

        random_signal = np.array([(v / len_chunk) / (mc / len_random_chunk) ** 2 for v in pc])
        random_signal = spectral.mean_smooth(random_signal, 3)

        # detrended smooth autocorrelation
        initial_values_dict = {'a_0': np.mean(random_signal[4:]), 'a_1': 0., 'a_2': 0.}
        params, obj_func = non_linear.create_quadratic_model(initial_values_dict)
        x = np.arange(len(random_signal[4:]))
        non_linear_fitter = non_linear.NonLinearFitting(obj_func, params, x, random_signal[4:])
        _, predicted = non_linear_fitter.least_squares()
        random_signal_detrended = random_signal[4:] - predicted

        # DTFT
        if counter == 1:
            label = 'randomized'
        else:
            label = None
        dtft_spectrum_plot(random_signal_detrended, ax_spec, title='periodogram',
                           norm=False, color='grey', alpha=0.5, label=label)

    dtft_spectrum_plot(signal_detrended, ax_spec, title='periodogram', norm=False, label='observed')
    ax_spec.legend()
    plt.rcParams['savefig.facecolor'] = fig_spec.get_facecolor()
