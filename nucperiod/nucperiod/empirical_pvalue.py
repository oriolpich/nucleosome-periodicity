#=================================================================================
# Name          : empirical_pvalue.py
# Author        : Oriol Pich
# Date          : 21/03/18
# Description   :
#=================================================================================

import numpy as np

from nucperiod.spectral import compute

def minor_in(mutation_array, expected_array, rand):

    difference = (np.array(mutation_array) - np.array(expected_array)) / np.array(expected_array)
    x_s, y_s, snr_muts, peak_muts = compute(difference, center=10.15, low_t=0,
                                                     high_t=115, low_p=5, high_p=20,
                                                     norm=True)
    count_snr_higher = 0

    # loop over each randomization
    for r in rand:

        # expected SNR
        difference = (np.array(r)-np.array(expected_array))/np.array(expected_array)

        x_s, y_s, snr_random, peak_random = compute(difference, center=10.15, low_t=0,
                                                             high_t=115, low_p=5, high_p=20,
                                                              norm=True)
        # get the empirical SNR pval
        if snr_random > snr_muts:
            count_snr_higher += 1

    pval_snr = count_snr_higher / len(rand)

    return  pval_snr


def zoomout(mutation_array, expected_array, rand):


    # observed SNR
    difference = (np.array(mutation_array) - np.array(expected_array)) / np.array(expected_array)

    x_s, y_s, snr_muts, peak_muts = compute(difference,center=191.3, low_t=0, high_t=1996, low_p=50,
                                                     high_p=250, norm=True)


    count_snr_higher = 0

    for r in rand:

        # expected SNR
        difference = (np.array(r) - np.array(expected_array)) / np.array(expected_array)

        x_s, y_s, snr_random, peak_random = compute(difference, center=191.3, low_t=0, high_t=1996,
                                                             low_p=50, high_p=250, norm=True)

        if snr_random > snr_muts:
            count_snr_higher += 1


    pval_snr = count_snr_higher / len(rand)

    return pval_snr
