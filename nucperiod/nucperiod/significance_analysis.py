#=================================================================================
# Name          : significance_analysis.py
# Author        : Oriol Pich
# Date          : 5/04/18
# Description   :
#=================================================================================


import numpy as np
from nucperiod.positions import ALL_MINORS_IN, ALL_MINORS_OUT


def minor_in(xvals, observed, expected):

    expected_minor_in = 0
    observed_minor_in = 0
    expected_minor_out = 0
    observed_minor_out = 0

    for ix, xpos in enumerate(xvals):
        if xpos in ALL_MINORS_IN:
            observed_minor_in += observed[ix]
            expected_minor_in += expected[ix]
        if xpos in ALL_MINORS_OUT:
            observed_minor_out += observed[ix]
            expected_minor_out += expected[ix]

    factor_norm = (observed_minor_in + observed_minor_out) / (expected_minor_in + expected_minor_out)
    expected_minor_in = expected_minor_in*factor_norm
    expected_minor_out = expected_minor_out*factor_norm

    return observed_minor_in, observed_minor_out, expected_minor_in, expected_minor_out


def zoomout(observed, expected):

    nucleosome_in_observed = np.sum(observed[970:1031])
    nucleosome_in_expected = np.sum(expected[970:1031])

    linker_in_observed_left = np.sum(observed[872:903])
    linker_in_expected_left = np.sum(expected[872:903])
    linker_in_observed_right = np.sum(observed[1093:1124])
    linker_in_expected_right = np.sum(expected[1093:1124])

    linker_observed = linker_in_observed_left + linker_in_observed_right
    linker_expected = linker_in_expected_left + linker_in_expected_right

    factor_norm = (nucleosome_in_observed + linker_observed) / (linker_expected + nucleosome_in_expected)

    expected_in_nucleosome = nucleosome_in_expected*factor_norm
    expected_in_linker = linker_expected*factor_norm

    return nucleosome_in_observed, linker_observed, expected_in_nucleosome, expected_in_linker
