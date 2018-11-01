from collections import defaultdict

import numpy as np
import pandas as pd

ALL_MINORS_IN = set()
ALL_MINORS_OUT = set()

# Get the positions of minor in and minor out zones
ROUNDED_OUTSIDE = [0, -10, -21, -31, -41, -51, -62,
                   10, 21, 31, 41, 51, 62]

ROUNDED_INSIDE = [-5, -15, -26, -36, -46, -57, -67,
                  5, 15, 26, 36, 46, 57, 67]

for r in ROUNDED_INSIDE:
    for i in range(r - 2, r + 3):
        ALL_MINORS_IN.add(i)

for r in ROUNDED_OUTSIDE:
    for i in range(r - 2, r + 3):
        ALL_MINORS_OUT.add(i)


ALL_MINORS_IN.add(-18)
ALL_MINORS_IN.add(18)
ALL_MINORS_OUT.add(-54)
ALL_MINORS_OUT.add(54)


def get_outphase(val):
    """
    Classify whether the reads are inphase or outphase
    """
    if val in ALL_MINORS_OUT:
        cl = 'INPHASE'
    elif val in ALL_MINORS_IN:
        cl = 'OUTPHASE'
    else:
        cl = np.nan

    return cl


def annotate_midpoints(df):
    # get the distance from the center
    df['distance'] = df['e1'] - df['pos2']

    # classify reads
    df['class'] = df['distance'].apply(get_outphase)

    annotated_dyads = defaultdict(list)

    # get the score for each nucleosome
    for i, data in df.groupby(by='ID'):
        # number of reads in phase (untested)
        inphase = data[data['class'] == 'INPHASE']['reads'].sum()
        total = data['reads'].sum()

        score = inphase / total

        chr_, pos1, pos2 = i.split('_')

        annotated_dyads['ID'].append(i)
        annotated_dyads['chr'].append(chr_)
        annotated_dyads['pos1'].append(int(pos1))
        annotated_dyads['pos2'].append(int(pos2))
        annotated_dyads['score_rot'].append(score)
        annotated_dyads['reads_in'].append(inphase)
        annotated_dyads['total_reads'].append(total)

    return pd.DataFrame(annotated_dyads)
