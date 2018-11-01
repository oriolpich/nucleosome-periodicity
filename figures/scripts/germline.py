import matplotlib.pyplot as plt
from nucperiod import plot as nucperiod_plt


def sapiens_bars(snr_high, snr_medium, snr_low):
    nucperiod_plt.config_params(7)

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(1.3, 1.5))
    xvals = [0, 1, 2]
    yvals = [snr_low, snr_medium, snr_high]

    axs.bar(xvals, yvals, label=['low', 'medium', 'high'], color=['#afdde9ff', '#afdde9ff', '#afdde9ff'])

    axs.set_xticks([0, 1, 2])
    axs.set_xticklabels(('very-low', 'low' ,'high'), rotation=90)
    axs.set_ylabel('SNR')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.tight_layout()
