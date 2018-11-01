import os

import click
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import ticker
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from matplotlib.ticker import FormatStrFormatter

from nucperiod import positions, plot, spectral, empirical_pvalue, significance_analysis

pandas2ri.activate()

# Configure plot parameters
plot.config_params_full()

# for y-axis formating
FORMATTER = ticker.ScalarFormatter(useMathText=True)
FORMATTER.set_scientific(True)
FORMATTER.set_powerlimits((-1, 1))


def first_observed_expected_plot(axs, xvals, norm_exp, yvals, half_win):

    spline = robjects.r["smooth.spline"]

    # plot expected
    axs.plot(xvals, norm_exp, linewidth=0.3, color='grey', alpha=0.5)

    # plot observed
    axs.plot(xvals, yvals, linewidth=0.3, color='red', alpha=0.5)

    # plot smoothed expected
    val = spline(xvals, norm_exp)
    axs.plot(list(val[0]), list(val[1]), label='Expected', color='#3d3d3d', linewidth=1.2, alpha=0.8)

    # plot smoothed observed
    val = spline(xvals, yvals)
    axs.plot(list(val[0]), list(val[1]), label='Observed', color='#e83e3e', linewidth=1.2, alpha=0.8)

    axs.set_ylabel('Mutation rate')
    axs.set_xlabel('Distance from dyad (bp)')

    axs.yaxis.set_major_formatter(FORMATTER)

    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    axs.spines['bottom'].set_visible(False)
    axs.get_xaxis().set_visible(False)

    axs.set_xlim(-half_win - 1, half_win + 1)


def plot_next_nucleosomes(good_local_maxima, axs, half_win, first_ax):

    # plot squares
    ax2 = axs.twinx()
    ax2.set_xticks(axs.get_xticks())

    for xd in good_local_maxima:
        axs.axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)
        first_ax.axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)

    for ix in range(len(good_local_maxima) - 1):
        maxima = good_local_maxima[ix]

        # this covers half of the nucleosome to the right
        start_next = maxima + 73

        # this gets the putative length of the linker
        length = (good_local_maxima[ix + 1] - 73) - start_next

        # plot the linker patch
        ax2.add_patch(plt.Rectangle((start_next, -0.105), length, 0.06, color='#7570b3',
                                    clip_on=False, linewidth=0, ))

    # plot the local maxima patches
    for maximum in good_local_maxima:
        ax2.add_patch(plt.Rectangle((maximum - 73, -0.105), 147, 0.06, color='#d95f02',
                                    clip_on=False, linewidth=0, ))

    # move axes down
    axs.spines["bottom"].set_position(("axes", -0.1093))

    ax2.set_xlim(-half_win - 1, half_win + 1)
    axs.set_xlim(-half_win - 1, half_win + 1)
    first_ax.set_xlim(-half_win - 1, half_win + 1)

    ax2.get_yaxis().set_ticks([])

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)


def second_smoothed_difference_plot(axs, xvals, difference, half_win, first_ax, dyads_closer_file):

    spline = robjects.r["smooth.spline"]

    # plot difference
    axs.plot(xvals, difference, linewidth=0.3, color='lightblue', alpha=0.7)
    axs.set_xlabel('Distance from dyad (bp)')
    axs.yaxis.set_major_formatter(FORMATTER)

    # plot smoothed expected
    try:
        val = spline(xvals, difference, df=60)
    except:
        pass
    else:

        if half_win <= 147:

            axs.set_ylabel('Relative increase')

            # paint the minors, etc
            plot.wave_painting_zoomin(list(val[0]), list(val[1]), axs, '#515e19ff', '#f3bb00', '#AFC6E9', 1.2)
            # plot the rectangles
            plot.rectangles_drawing(axs, half_win, '#f3bb00', '#515e19ff')

        else:
            axs.set_ylabel('Relative increase')

            # decorate the plots with next nucleosomes if large window
            good_local_maxima = np.load(dyads_closer_file)

            # pint the nucleosomes and dyads, etc
            plot.wave_painting_zoomout(list(val[0]), list(val[1]), axs, good_local_maxima)

            # plot the rectangles
            plot_next_nucleosomes(good_local_maxima, axs, half_win, first_ax)
            axs.set_ylim(-0.2, 0.2)

    # remove math annotation
    axs.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=False))
    axs.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)


def period_spectra(axs, difference, zoomin):

    if zoomin:
        # non normalized to get the maximum amplitude
        x, y, snr, peak = spectral.compute(difference, low_t=0, high_t=115, low_p=5, high_p=20, norm=True)

        axs.plot(x, y, color='#31a354', lw=1.2)
        axs.set_xlim(4, 20)

    else:  # if the window is big, we change the scope of the periodicity

        # non normalized to get the maximum amplitude
        x, y, snr, peak = spectral.compute(difference, low_t=0, high_t=1996, low_p=50, high_p=250, norm=True)

        # plot the spectogram
        axs.plot(x, y, color='#31a354', lw=1.2)
        axs.set_ylim(0, 10)
        axs.set_xlim(100, 250)

    # this will be the header
    out_t = ' '
    axs.set_title(out_t, fontsize=7.5)
    axs.set_ylabel('Power (a.u.)')
    axs.set_xlabel('Period (bp)')

    return round(snr, 2), round(peak, 2)


def error_bars_expected(rand, norm_exp, total_nucleosomes, total_samples):
    """
    We will get the error bar per position defined by the standard deviation
    """
    # lower boundary
    lower_b = np.array(norm_exp) - (np.std(rand, axis=0)/total_nucleosomes) / total_samples

    # upper boundary
    higher_b = np.array(norm_exp) + (np.std(rand, axis=0)/total_nucleosomes) / total_samples

    return lower_b, higher_b


def load_obs_exp(mutations_file, dyads_closer_file, expected_nucleosomes_file, output_folder, zoomin, total_samples, total_nucleosomes):

    # determine the size of the windows
    if zoomin:
        half_win = 58
    else:
        half_win = 1000

    xvals = np.arange(-half_win, half_win + 1)

    # load the data we generated with the intersection of the mutations and nucleosomes
    df_mutations = pd.read_csv(mutations_file, sep='\t', names=['id', 'pos', 'kmer'])
    d_pos = df_mutations['pos'].value_counts().to_dict()

    # here we will add the observed mutation count
    real_observed = np.array([d_pos.get(i, 0) for i in range(-half_win, half_win + 1)])

    # normalize all counts
    total_muts = sum(real_observed)
    observed_y = real_observed / total_nucleosomes
    normalized_observed = observed_y / total_samples

    # load all the expected mutations
    expected_nucleosomes_data = np.load(expected_nucleosomes_file)
    exp_array = expected_nucleosomes_data['freq']
    normalized_nucl_exp_array = np.array(exp_array) / total_nucleosomes
    norm_exp = np.array(normalized_nucl_exp_array) / total_samples

    # load randomizations
    rand = []
    if 'rand' in expected_nucleosomes_data:
        rand = expected_nucleosomes_data['rand']

    nrands = len(rand)

    # create the plot structure
    if zoomin:
        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(4.88 / 2.13, 7.5 / 1.83))
    else:
        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(4.88 / 1.96, 7.5 / 1.79))

    plt.subplots_adjust(hspace=0.001)

    # write dashed vertical lines
    if zoomin:
        # plot in the first and second plots
        for xd in positions.DYAD_X_SMALL:
            axs[0].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)
            axs[1].axvline(xd, color='black', linestyle='--', lw=0.2, alpha=0.5)

    # RANDOMIZATION if any #
    if nrands > 0:
        # get empirical pvalues
        if zoomin:
            empirical_pval_snr = \
                empirical_pvalue.minor_in(real_observed, exp_array, rand)
        else:  # get empirical pvalues in zoomout
            empirical_pval_snr = \
                empirical_pvalue.zoomout(real_observed, exp_array, rand)

        # error bars of the distribution
        lower_b, higher_b = error_bars_expected(rand, norm_exp, total_nucleosomes, total_samples)

        # fill grey area with higher and lower boundaries which are equivalent of SDM
        axs[0].fill_between(xvals, lower_b, higher_b, color='#bdbdbd', alpha=0.55, linewidth=0.3)

    # for quick tests
    else:
        empirical_pval_snr = 'nan'

    ## FIRST PLOT #
    first_observed_expected_plot(axs[0], xvals, norm_exp, normalized_observed, half_win)

    ## SECOND PLOT ##
    difference = (normalized_observed - norm_exp) / norm_exp

    second_smoothed_difference_plot(axs[1], xvals, difference, half_win, axs[0], dyads_closer_file)

    # gtest analysis
    if zoomin:
        observed_minor_in, observed_minor_out, expected_minor_in, expected_minor_out = \
            significance_analysis.minor_in(xvals, real_observed, exp_array)

    # special gtest nucleosome-linker
    else:
        observed_minor_in, observed_minor_out, expected_minor_in, expected_minor_out = \
            significance_analysis.zoomout(real_observed, exp_array)


    ## THIRD PLOT (SPECTROGRAM)
    snr, peak = period_spectra(axs[2], difference, zoomin)

    # crosscorrelation analysis
    N = half_win * 2 + 1
    if zoomin:
        P = 10.3
        bounds = [-58, 59]
    else:
        P = 191.3
        bounds = [-1000, 1001]

    current_phase = spectral.assess_phase(xvals, difference, N, P, bounds, zoomin)
    if empirical_pval_snr == 0:
        empirical_pval_snr_name = ' < 0.001'
    else:
        empirical_pval_snr_name = ' = {}'.format(empirical_pval_snr)
    axs[1].set_title('n = {}, SNR = {}\nMP = {}, p-value{}, phase = {}'.format(total_samples, snr, peak,
                                                                               empirical_pval_snr_name, current_phase), fontsize=7.5)

    axs[2].spines['right'].set_visible(False)
    axs[2].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['left'].set_visible(False)

    ## write all analyses done
    outf = os.path.join(output_folder, 'obsexp.tsv')

    with open(outf, 'wt') as outfile:
        out = '\t'.join(map(str, [total_muts, observed_minor_in, observed_minor_out,
                         expected_minor_in, expected_minor_out, snr, peak,
                         current_phase, empirical_pval_snr]))

        outfile.write(out)
        outfile.write('\n')

    plt.tight_layout()

    outf = os.path.join(output_folder, 'obsexp.png')
    plt.savefig(outf, dpi=600, clip=False, bbox_inches='tight')
    outf = os.path.join(output_folder, 'obsexp.svg')
    plt.savefig(outf)

    plt.close()


@click.command()
@click.argument('mutation_file', metavar='<MUTATIONS>')
@click.argument('expected_nucleosomes_file', metavar='<EXPECTED NUCLEOSOMES>')
@click.argument('output', metavar='<OUT FOLDER>')
@click.option('--zoomin/--zoomout', help='Zoom in and zoom out analysis')
@click.option('--dyads', 'dyads_closer_file', type=click.Path(), default=None, help='Dyads closer file (numpy list). Required for zoomout')
@click.option('--samples', type=click.INT, help='Number of samples')
@click.option('--nucleosomes', type=click.INT, help='Number of nucleosomes')
def cli(mutation_file, dyads_closer_file, expected_nucleosomes_file, output, zoomin, samples, nucleosomes):
    """

    :param mutation_file:
    :param dyads_closer_file:
    :param expected_nucleosomes_file:
    :param output:
    :param zoomin:
    :param samples:
    :param nucleosomes:
    :return:
    """

    if not zoomin and dyads_closer_file is None:
        raise click.BadParameter('Dyad file is required for zoomout')

    load_obs_exp(mutation_file, dyads_closer_file, expected_nucleosomes_file, output, zoomin, samples, nucleosomes)


if __name__ == '__main__':
    cli()
