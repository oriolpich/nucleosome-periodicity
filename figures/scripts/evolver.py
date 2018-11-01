import gzip
import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
from nucperiod import plot as nucperiod_plt, positions, spectral
from rpy2 import robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()


FORMATTER = ticker.ScalarFormatter(useMathText=True)
FORMATTER.set_scientific(True)
FORMATTER.set_powerlimits((-3, 3))


def reconstruct_seq_centered(seq, nucleosome_pos):
    """ Reconstruct the sequence and get the peak and SNR for each iteration """

    # equivalence
    d_nucleotide = {0: 1,  # 'A',
                    1: 0,  # 'C',
                    2: 1,  # 'T',
                    3: 0,  # 'G'
                    }

    seqd = np.vectorize(d_nucleotide.get)(seq)
    array_nuc = []

    # select only the nucleosome positions
    for pos in nucleosome_pos:
        array_nuc.append(seqd[pos - 58:pos + 59])

    # do the stack
    center_nucleosome = np.sum(array_nuc, axis=0) / len(array_nuc)

    return center_nucleosome


def sequence(evolver):
    nucperiod_plt.config_params_full(7)

    ww = reconstruct_seq_centered(evolver[1], evolver[2])

    spline = robjects.r["smooth.spline"]

    xvals = np.arange(-58, 59)
    div = ww
    yvals = ww
    val = spline(np.arange(-58, 59), ww, df=60)

    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(4.88 / 2.13, 7.5 / 1.83))
    plt.subplots_adjust(hspace=0.001)

    for xd in positions.DYAD_X_SMALL:
        ax[0].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)
        ax[1].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)

    ax[0].set_title('Sequence', fontsize=7.5)
    ax[0].plot(xvals, yvals, c='#1f78b4ff', linewidth=0.8, )

    ax[0].yaxis.set_major_formatter(FORMATTER)

    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[0].get_xaxis().set_visible(False)
    ax[0].set_xlim(-58, 59)
    ax[0].set_ylabel('Damage reads')

    ax[1].plot(xvals, div, linewidth=0.3, color='lightblue', alpha=0.7)

    nucperiod_plt.wave_painting_zoomin(list(val[0]), list(val[1]), ax[1], color_out='#add8e6ff', color_in='#8b0000ff',
                                       color_other='#add8e6ff', line_width=1.2)
    nucperiod_plt.rectangles_drawing(ax[1], 58, '#ffcf12', '#728130')

    x, y, snr, peak = spectral.compute(list(val[1]), center=10.3, low_t=0, high_t=115, low_p=5, high_p=20, norm=True)

    ax[2].plot(x, y, color='#31a354', lw=1.2)
    ax[2].set_xlim(4, 20)

    ax[1].set_title('n = 1, SNR = {}\nMP = {}, phase = 1'.format(round(snr, 2), round(peak, 2),
                                                                 ), fontsize=7.5)
    ax[2].set_title('SNR = {}, MP = {}, phase = 1'.format(round(snr, 2), round(peak, 2)), fontsize=7.5)
    ax[1].set_ylabel('WW proportion')
    ax[1].set_xlabel('Distance from dyad (bp)')
    ax[2].set_ylabel('Power')
    ax[2].set_xlabel('Period (bp)')
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['left'].set_visible(False)

    plt.tight_layout()


def iterations(evolver1, evolver2, evolver3):
    nucperiod_plt.config_params_full(7)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3,2.7))
    plt.plot(evolver1[3], lw=2, c='#ffa500ff')
    plt.plot(evolver2[3], c='#e3000064', alpha=0.3)
    plt.plot(evolver3[3], marker='o', markersize=1, c='#0000ffff', alpha=0.3)

    plt.xlim(0, 60)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks(np.arange(0, 70, 10), [i*100 for i in np.arange(0, 70, 10)])
    plt.ylabel('SNR 10bp WW periodicity')
    plt.xlabel('Iterations (time)')


def iterations_all(evolver1_files, evolver2_files, evolver3_files):
    nucperiod_plt.config_params_full(7)

    real = []
    noprob = []
    nobias = []

    for f1, f2, f3 in zip(evolver1_files, evolver2_files, evolver3_files):
        t = pickle.load(gzip.open(f1))
        real.append(t[3])
        t = pickle.load(gzip.open(f2))
        noprob.append(t[3])
        t = pickle.load(gzip.open(f3))
        nobias.append(t[3])

    mean_all = np.mean(np.array(real).T, axis=1)
    mean_noprob = np.mean(np.array(noprob).T, axis=1)
    mean_nobias = np.mean(np.array(nobias).T, axis=1)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 2.7))

    plt.plot(np.arange(len(t[3])), mean_all, lw=4, c='#ffa500ff')

    upper = []
    lower = []
    for pos in np.array(real).T:
        upper.append(np.percentile(pos, 90))
        lower.append(np.percentile(pos, 10))

    plt.fill_between(np.arange(len(t[3])), lower, upper, color='grey', alpha=0.1)

    plt.plot(np.arange(len(t[3])), mean_noprob, c='#e3000064')

    upper = []
    lower = []
    for pos in np.array(noprob).T:
        upper.append(np.percentile(pos, 90))
        lower.append(np.percentile(pos, 10))

    plt.fill_between(np.arange(len(t[3])), lower, upper, color='grey', alpha=0.1)

    plt.plot(np.arange(len(t[3])), mean_nobias, '--', c='#0000ffff')

    upper = []
    lower = []
    for pos in np.array(nobias).T:
        upper.append(np.percentile(pos, 90))
        lower.append(np.percentile(pos, 10))

    plt.fill_between(np.arange(len(t[3])), lower, upper, color='grey', alpha=0.1)

    plt.xlim(0, 60)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks(np.arange(0, 70, 10), [i*100 for i in np.arange(0, 70, 10)])
    plt.ylabel('SNR 10bp WW periodicity')
    plt.xlabel('Iterations (time)')
