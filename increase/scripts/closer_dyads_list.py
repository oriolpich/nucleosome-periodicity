
import logging

import click
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from scipy.signal import argrelextrema

pandas2ri.activate()


def get_next_dyads_locations(dyads, outfile):

    # read the intersected file
    df = pd.read_csv(dyads, sep='\t', names=['ID', 'distance'])

    # get the distance excluding the nucleosomes that mapp themselves
    df = df[(df['distance']!=0)]

    # create the dictionary
    dic_distance = df['distance'].value_counts().to_dict()

    # sort the dictionary and get info of all points
    half_win = 1000
    yvals_closer_dyads= []
    xvals_closer_dyads = []
    for ix, i in enumerate(range(-half_win, half_win + 1)):
        yvals_closer_dyads.append(dic_distance.get(i, 0))
        xvals_closer_dyads.append(i)

    # smooth the count
    spline = robjects.r["smooth.spline"]
    normalized_dyad_count = np.array(yvals_closer_dyads) / np.sum(yvals_closer_dyads)
    val = spline(xvals_closer_dyads, normalized_dyad_count, df = 80)

    # find local maxima
    greater_points = np.array(val[0])[argrelextrema(np.array(list(val[1])), np.greater)[0]]
    good_local_maxima = [i for i in greater_points if abs(i) > 147 and abs(i) < 900]

    # add zero
    good_local_maxima.append(0)

    # sort the list
    good_local_maxima = sorted(good_local_maxima)
    np.save(outfile, good_local_maxima)


@click.command()
@click.argument('dyads', metavar='DYADS', type=click.Path(exists=True))
@click.argument('outfile', metavar='OUTFILE', type=click.Path())
def cli(dyads, outfile):
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG, datefmt='%H:%M:%S')
    get_next_dyads_locations(dyads, outfile)


if __name__ == '__main__':
    cli()
