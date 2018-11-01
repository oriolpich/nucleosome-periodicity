import collections
import gzip

import matplotlib.pyplot as plt
import numpy as np
from bgreference import refseq
from rpy2 import robjects
from rpy2.robjects import pandas2ri

from nucperiod import plot as nucperiod_plot, positions, spectral

pandas2ri.activate()


def proportion_ww(genome, file):
    wanted = 'AT'

    count_seq = 0
    data = collections.defaultdict(int)

    with gzip.open(file, 'rt') as infile:
        next(infile)
        for line in infile:
            line_spl = line.rstrip().split('\t')

            try:
                seq = refseq(genome, line_spl[0], int(line_spl[2]) - 73, 148)
                count_seq += 1

                for i in range(len(seq)):
                    if (seq[i] in wanted) & (seq[i + 1] in wanted):
                        data[i] += 1
            except:
                continue

    return data, count_seq


def sequence(ww_proportion, total):

    vals = []
    xvals = []
    for ix, v in enumerate(range(-73, 74)):
        vals.append(ww_proportion[ix])
        xvals.append(v)
    vals_s = vals[15:132]
    return np.array(vals_s) / total


def plot(groups):
    nucperiod_plot.config_params_full()
    # Compute proportions
    proportion_normalized = collections.defaultdict(lambda: collections.defaultdict(dict))

    for species in groups.keys():
        genome, file = groups[species]
        ww, total = proportion_ww(genome, file)
        proportion_normalized[species] = sequence(ww, total)

    # Do the plot
    color_out = '#add8e6ff'
    color_in = '#8b0000ff'
    half_win = 58

    fig, ax = plt.subplots(nrows=1, ncols=len(groups.keys()), figsize=(16, 1.8))

    xvals = [i for i in range(-58, 59)]

    for ix, species in enumerate(groups.keys()):
        ax[ix].set_ylabel('WW proportion')

        for xd in positions.DYAD_X:
            ax[ix].axvline(xd, color='black', linestyle='--', lw=0.5, alpha=0.5)
        spline = robjects.r["smooth.spline"]
        yvals = proportion_normalized[species]
        val = spline(xvals, list(yvals))

        x, y, snr, peak = spectral.compute(list(val[1]), low_t=0, high_t=115, low_p=5, high_p=20, norm=True)
        title = '$\it{}$\n\nSNR = {}, MP = {}'.format(species, round(snr, 3), round(peak, 3))
        ax[ix].set_title(title)

        nucperiod_plot.wave_painting_zoomin(list(val[0]), list(val[1]), ax[ix], color_in, color_out, color_out, 1.6)
        nucperiod_plot.rectangles_drawing(ax[ix], half_win, '#ffcf12', '#728130')


    plt.tight_layout()
