import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import ticker
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

from nucperiod import plot as nucperiod_plot, positions, spectral

# for y-axis formating
FORMATTER = ticker.ScalarFormatter(useMathText=True)
FORMATTER.set_scientific(True)
FORMATTER.set_powerlimits((-3, 3))


spline = robjects.r["smooth.spline"]


def plot(data):
    nucperiod_plot.config_params_full()

    df_polarized = pd.read_csv(data, sep='\t', names=['ID', 'pos', 'ref', 'alt'])

    CT_muts = df_polarized[((df_polarized['ref'] == 'C') & (df_polarized['alt'] == 'T')) | ((df_polarized['ref'] == 'G') & (df_polarized['alt'] == 'A'))]
    df_polarized = df_polarized[((df_polarized['ref'] == 'C') | (df_polarized['ref'] == 'G'))]

    d_muts = CT_muts['pos'].value_counts().to_dict()
    d_pols = df_polarized['pos'].value_counts().to_dict()

    muts = []
    pol = []
    xvals = []

    for i in range(-58, 59):
        muts.append(d_muts[i])
        pol.append(d_pols[i])
        xvals.append(i)

    div = np.array(muts) / np.array(pol)

    val = spline(xvals, div)

    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(4.88 / 2.13, 7.5 / 1.83))
    plt.subplots_adjust(hspace=0.001)

    for xd in positions.DYAD_X_SMALL:
        ax[0].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)
        ax[1].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)

    ax[0].set_title('Ancestral', fontsize = 7.5)

    ax[0].plot(xvals, div, linewidth=0.3, color='lightblue', alpha=0.7)

    nucperiod_plot.wave_painting_zoomin(list(val[0]), list(val[1]), ax[1], '#515e19ff', '#f3bb00', '#AFC6E9', 1.2)

    ax[1].plot(xvals, div, linewidth=0.3, color='lightblue', alpha=0.7)

    nucperiod_plot.wave_painting_zoomin(list(val[0]), list(val[1]), ax[1], '#515e19ff', '#f3bb00', '#AFC6E9', 1.2)
    nucperiod_plot.rectangles_drawing(ax[1], 58, '#f3bb00', '#515e19ff')

    x, y, snr, peak = spectral.compute(list(val[1]), center=10.3,
                                       low_t=0, high_t=115, low_p=5, high_p=20, norm=True)

    ax[2].plot(x, y, color='#31a354', lw=1.2)
    ax[2].set_xlim(4, 20)

    ax[1].set_title('n = 1, SNR = {}\nMP = {}, p-value = n1n, phase = -1'.format(round(snr, 2), round(peak, 2),
                                                                                 ), fontsize=7.5)
    ax[2].set_title('SNR = {}, MP = {}, phase = -1'.format(round(snr, 2), round(peak, 2)), fontsize=7.5)
    ax[1].set_ylabel('Divergence')
    ax[1].set_xlabel('Distance from dyad (bp)')
    ax[2].set_ylabel('Power')
    ax[2].set_xlabel('Period (bp)')
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['left'].set_visible(False)

    plt.tight_layout()
