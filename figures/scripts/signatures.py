import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

from nucperiod import plot as nuceriod_plt, increase


np.random.seed(123)


COLORS = {
        'Signature_5': '#f9e612',
        'Signature_15': '#b15928',
        'Signature_2': '#ffff99',
        'Signature_4': '#6a3d9a',
        'Signature_9': '#cab2d6',
        'Signature_10': '#e31a1c',
        'Signature_18': '#c69f04',
        'Signature_26': '#fdbf6f',
        'Signature_17': '#33a02c',
        'Signature_7': 'black',
        'Signature_13': '#fb9a99',
        'Signature_28': '#b2df8a',
        'Signature_1': '#1f78b4',
        'Signature_27': '#a6cee3',
        'Signature_16': '#aaaa00',
        'Signature_14': '#aa00ff',
        'Signature_3': '#38aa9d',
        'Signature_6': '#9992aa',
        'Signature_12': '#aaa1a1',
        'Signature_30': '#7d3a3b',
        'Signature_11': 'green',
        'Signature_19': 'grey',
        'Signature_20': 'pink',
        'Signature_21': 'blue',
        'Signature_22': 'white',
        'Signature_23': 'darkblue',
        'Signature_24': 'orange',
        'Signature_25': 'darkorange',
        'Signature_29': 'grey',
        'Signature_8': '#fff7bc'
    }


def _load(files):
    df = increase.load_d(files)

    df['increase_in'] = df['observed_in'] - df['expected_in']
    df['prop_increase_in'] = df['increase_in'] / df['expected_in']
    df['outname'] = df['name'].str.replace('Signature_', 'Sign ')

    df['empirical_pvalue_snr'].replace(0, 0.001, inplace=True)
    df['empirical_pvalue_snr'] = df['empirical_pvalue_snr'].fillna(1)

    qvals = multipletests(df['empirical_pvalue_snr'].tolist(), method='fdr_bh')
    df['qvals_snr'] = qvals[1]

    return df


def zoomin(files):
    nuceriod_plt.config_params(11)

    df = _load(files)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 7), sharex=True)

    df.sort_values(by='qvals_snr', ascending=False, inplace=True)

    for i, row in df.iterrows():
        edgecolor = None
        snr = row['snr']
        if row['qvals_snr'] < 0.1 and snr > 8:
            edgecolor = 'darkred'

        if row['cross_validation_max'] < 0:

            ax[1].scatter(row['peak'], -np.log2(snr), s=200, linewidth=2, edgecolor=edgecolor, color='white',
                          label=row['name'], alpha=1)

            ax[1].scatter(row['peak'], -np.log2(snr), s=80, edgecolor='grey', linewidth=0.5, color=COLORS[row['name']],
                          label=row['name'], alpha=1)
            if edgecolor == 'darkred':
                ax[1].text(row['peak'] + 1, -np.log2(snr) - 0.15, row['outname'])

    for i, row in df.iterrows():
        edgecolor = None
        snr = row['snr']
        if row['qvals_snr'] < 0.1 and snr > 8:
            edgecolor = 'darkred'
        if row['cross_validation_max'] > 0:

            ax[0].scatter(row['peak'], np.log2(snr), s=200, linewidth=2, edgecolor=edgecolor, color='white',
                          label=row['name'], alpha=1)

            ax[0].scatter(row['peak'], np.log2(snr), s=80, edgecolor='grey', linewidth=0.5, color=COLORS[row['name']],
                          label=row['name'], alpha=1)

            if edgecolor == 'darkred':
                ax[0].text(row['peak'] + 1, np.log2(snr) - 0.15, row['outname'])

    plt.xlabel('Period')

    xlim = [i for i in range(8, 22, 2)]
    ax[0].set_xticks(xlim)
    ax[1].set_xticks(xlim)
    yvals = [i for i in range(2, 10, 2)]
    ax[0].set_yticks(yvals)

    yvals = [i for i in range(-8, 0, 2)]
    ax[1].set_yticks(yvals)
    ylabels = [str(2 ** abs(i)) for i in range(2, 10, 2)]
    ax[0].set_yticklabels(ylabels)

    ylabels = ['{}'.format(str(2 ** abs(i))) for i in range(-8, 0, 2)]
    ax[1].set_yticklabels(ylabels)

    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    ax[1].xaxis.set_ticks_position('top')
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['right'].set_visible(False)

    ax[0].set_ylim(1.5, 10)
    ax[1].set_ylim(-10, -1.5)

    plt.tight_layout()


def zoomout(files):
    nuceriod_plt.config_params(14)

    df = _load(files)

    df = df[df['nmuts_whole_nucleosome'] > 10000]

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 7), sharex=True)

    for i, row in df.iterrows():
        edgecolor = None
        snr = row['snr']
        if row['qvals_snr'] < 0.05 and snr > 8:
            edgecolor = 'darkred'

        if row['cross_validation_max'] > 0:

            ax[0].scatter(row['peak'], np.log2(snr), s=200, linewidth=2, edgecolor=edgecolor, color='white',
                          label=row['name'], alpha=1)

            ax[0].scatter(row['peak'], np.log2(snr), s=80, edgecolor='grey', linewidth=0.5, color=COLORS[row['name']] ,
                          label=row['name'], alpha=1)
            if edgecolor == 'darkred':
                ax[0].text(row['peak'] + 10, np.log2(snr) - 0.15, row['outname'])

    for i, row in df.iterrows():
        edgecolor = None
        snr = row['snr']
        if row['qvals_snr'] < 0.05 and snr>8:
            edgecolor = 'darkred'

        if row['cross_validation_max'] < 0:

            ax[1].scatter(row['peak'], -np.log2(snr), s=200, linewidth=2, edgecolor=edgecolor, color='white',
                          label=row['name'], alpha=1)

            ax[1].scatter(row['peak'], -np.log2(snr), s=80, edgecolor='grey', linewidth=0.5, color=COLORS[row['name']],
                          label=row['name'], alpha=1)

            if edgecolor == 'darkred':
                ax[1].text(row['peak'] + 10, -np.log2(snr) - 0.15, row['outname'])

    plt.ylabel('log2(SNR)')
    plt.xlabel('Period')

    xlim = [i for i in range(50, 260, 20)]
    ax[0].set_xticks(xlim)
    ax[1].set_xticks(xlim)
    yvals = [i for i in range(2, 10, 2)]
    ax[0].set_yticks(yvals)

    yvals = [i for i in range(-8, 0, 2)]
    ax[1].set_yticks(yvals)
    ylabels = [str(2 ** abs(i)) for i in range(2, 10, 2)]
    ax[0].set_yticklabels(ylabels)

    ylabels = [str(2 ** abs(i)) for i in range(-8, 0, 2)]
    ax[1].set_yticklabels(ylabels)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    ax[1].xaxis.set_ticks_position('top')
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['right'].set_visible(False)

    ax[0].set_ylim(1.5, 10)
    ax[1].set_ylim(-10, -1.5)

    plt.tight_layout()


def generate_table(files):
    df = _load(files)
    table_out = df[['name', 'peak', 'snr', 'empirical_pvalue_snr', 'qvals_snr', 'cross_validation_max']].copy()
    table_out.columns = ['Signature', 'Peak', 'SNR', 'P-value', 'Q-value', 'Phase']
    table_out.sort_values(by='Signature', inplace=True)
    return table_out


def rotational_bars(rot_high_files, rot_low_files):
    nuceriod_plt.config_params(11)

    df_high = increase.load_d(rot_high_files)
    df_high['rot'] = 'high'

    df_low = increase.load_d(rot_low_files)
    df_low['rot'] = 'low'

    df = pd.concat([df_high, df_low])

    fig, axs = plt.subplots(nrows=len(df.groupby(by='name')), ncols=1, figsize=(1.75, 8))
    order = ['low', 'high']

    for ix, (sig, data) in enumerate(df.groupby(by='name')):
        xvals = []
        yvals = []
        colors = []
        count = 0
        for i in order:
            val = data[data['rot'] == i]['snr'].tolist()[0]
            yvals.append(val)
            xvals.append(count)
            colors.append(COLORS[sig])
            count += 1
        axs[ix].bar(xvals, yvals, color=colors, label=['low', 'high'])
        axs[ix].set_xticks([0, 1])
        axs[ix].set_xticklabels(('low', 'high'), fontsize=11)
        axs[ix].set_ylabel('SNR')
        axs[ix].spines['right'].set_visible(False)
        axs[ix].spines['top'].set_visible(False)
    plt.tight_layout()


def compare(files_deconstructsigs, files_sigfit):
    nuceriod_plt.config_params(11)

    df_deconstructsigs = increase.load_d(files_deconstructsigs)
    df_deconstructsigs['control'] = 'deconstructsigs'

    df_sigfit = increase.load_d(files_sigfit)
    df_sigfit['control'] = 'sigfit'

    toplot = pd.concat([df_deconstructsigs, df_sigfit])

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 2))
    sig_list = []
    decon_list = []
    sig_list2 = []
    decon_list2 = []

    for i, data in toplot.groupby(by='name'):
        if len(data) == 2:
            row = data.iloc[0]
            snr_decon = data[data['control'] == 'sigfit']['snr'].iloc[0]
            snr_sigfit = data[data['control'] == 'deconstructsigs']['snr'].iloc[0]
            sig_list.append(snr_sigfit)
            decon_list.append(snr_decon)
            if row['cross_validation_max'] < 0:
                ax.scatter(-np.log(snr_decon), -np.log(snr_sigfit), c=COLORS[i])
                sig_list2.append(-np.log(snr_sigfit))
                decon_list2.append(-np.log(snr_decon))
            else:
                ax.scatter(np.log(snr_decon), np.log(snr_sigfit), c=COLORS[i])
                sig_list2.append(np.log(snr_sigfit))
                decon_list2.append(np.log(snr_decon))
    plt.xlabel('Period')

    ylabels = [str(2 ** abs(i)) for i in range(2, 10, 2)]
    yfinal = ylabels[::-1] + ylabels[1:]
    plt.xticks(np.arange(-6, 8, 2), yfinal)

    ylabels = [str(2 ** abs(i)) for i in range(2, 10, 2)]
    yfinal = ylabels[::-1] + ylabels[1:]
    plt.yticks(np.arange(-6, 8, 2), yfinal)

    slope, intercept, r_value, p_value, std_err = stats.linregress(sig_list2, decon_list2)
    xvals = np.arange(-6, 8, 2)
    yvals = [slope * y + intercept for y in xvals]
    plt.plot(xvals, yvals)
    R, pval = stats.pearsonr(sig_list, decon_list)

    plt.text(-4, 3, 'R = {}\npval = {}'.format(round(R, 3), round(pval, 3)))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel('SNR (Sigfit)')
    plt.xlabel('SNR (Deconstructsig)')

    plt.tight_layout()


def jitter(x):
    x = x + np.random.randn(1)*0.1
    return x


def cohorts(files, colors):

    nuceriod_plt.config_params_full(font_size=7)

    df = increase.load_d(files)

    df['increase_in'] = df['observed_in'] - df['expected_in']
    df['prop_increase_in'] = df['increase_in'] / df['expected_in']

    df['signature_name'] = df['name'].str.split('__').str.get(1)
    df['ttype'] = df['name'].str.split('__').str.get(0)

    df['empirical_pvalue_snr'].replace(0, 0.001, inplace=True)
    df['empirical_pvalue_snr'] = df['empirical_pvalue_snr'].fillna(1)

    qvals = multipletests(df['empirical_pvalue_snr'].tolist(), method='fdr_bh')
    df['qvals_snr'] = qvals[1]

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(7, 2))

    signature_groups = df.groupby('signature_name')

    for i in range(1, 31):
        ix = i -1
        try:
            data = signature_groups.get_group('Signature_{}'.format(i))
        except KeyError:
            continue
        else:

            for _, row in data.iterrows():

                if row['cross_validation_max'] > 0:
                    if row['qvals_snr'] < 0.05 and row['snr'] > 8 and 10 < row['peak'] < 10.4:
                        axs[0].scatter(jitter(ix), np.log2(row['snr']), c=colors[row['ttype']], s=25,
                                       edgecolor='black', linewidth=0.5)
                    else:
                        axs[0].scatter(jitter(ix), np.log2(row['snr']), c='grey', s=6)

                elif row['cross_validation_max'] < 0:

                    if row['qvals_snr'] < 0.05 and row['snr'] > 8 and 10 < row['peak'] < 10.4:
                        axs[1].scatter(jitter(ix), -np.log2(row['snr']), c=colors[row['ttype']], s=25,
                                       edgecolor='black', linewidth=0.5)
                    else:
                        axs[1].scatter(jitter(ix), -np.log2(row['snr']), c='grey', s=8)


    yvals = [i for i in range(2, 10, 2)]
    axs[0].set_yticks(yvals)

    yvals = [i for i in range(-8, 0, 2)]
    axs[1].set_yticks(yvals)
    ylabels = [str(2 ** abs(i)) for i in range(2, 10, 2)]
    axs[0].set_yticklabels(ylabels)

    ylabels = ['{}'.format(str(2 ** abs(i))) for i in range(-8, 0, 2)]
    axs[1].set_yticklabels(ylabels)

    axs[0].spines['right'].set_visible(False)
    axs[0].spines['top'].set_visible(False)
    axs[0].set_xlim(-1, 30)
    xlabels = [i for i in range(1, 31)]
    xpos = [i for i in range(0, 30)]

    axs[1].set_xticklabels(xlabels)
    plt.xticks(xpos, xlabels)
    axs[1].xaxis.set_ticks_position('top')
    axs[1].spines['bottom'].set_visible(False)
    axs[1].spines['right'].set_visible(False)

    axs[0].set_ylim(1.5, 10)
    axs[1].set_ylim(-10, -1.5)
    plt.tight_layout()
