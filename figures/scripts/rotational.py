import matplotlib.pyplot as plt

from nucperiod import plot as nucperiod_plt


def plot_bars(snr_high, snr_low):
    nucperiod_plt.config_params()

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(0.9,1.1,))

    xvals = [0, 1]
    yvals = [snr_low, snr_high]
    axs.bar(xvals, yvals, label=['low','high'], color=['#afdde9ff', '#afdde9ff'])

    plt.xticks(xvals, ['low', 'high'], fontsize=7)
    axs.set_xlabel('SNR', fontsize = 7)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.tight_layout()
