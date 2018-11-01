import collections
import functools
import json
import os
from multiprocessing import Pool

import click
import numpy as np
import pandas as pd

from nucperiod import spectral


def counts_to_spectrum(counts_json, low_t, high_t, low_p, high_p, resolution):
    with open(counts_json, 'rt') as f:
        counts_dict = json.load(f)
    name = os.path.basename(counts_json).replace('.json', '')
    pair_count = np.array(counts_dict['pair_count'])
    motif_count = counts_dict['motif_count']
    len_chunk = counts_dict['chunk_len']
    signal = np.array([(v / len_chunk) / (motif_count / len_chunk) ** 2 for v in pair_count])
    x = np.linspace(low_p, high_p, resolution)
    y = spectral.spectrum(signal, x, norm=False, low_t=low_t, high_t=high_t, low_p=low_p, high_p=high_p)
    return name, y


@click.command(context_settings={'help_option_names': ['h', '--help']})
@click.argument('input', type=click.Path(), metavar='<FOLDER>')
@click.argument('output', type=click.Path(), metavar='<FILE>')
@click.option('--cores', type=click.INT, help='cores to use', default=1)
@click.option('--randomized', is_flag=True, default=False, help='If the data comes from observed or randomized')
@click.option('--table', type=click.Path(), help='File to output results in tabular form', default=None)
def cli(input, output, cores, randomized, table):
    low_t, high_t = 30, 100
    low_p, high_p = 5, 20
    resolution = 100

    f = functools.partial(counts_to_spectrum, low_t=low_t, high_t=high_t, low_p=low_p, high_p=high_p, resolution=resolution)

    if randomized:  # search in nested directories
        spectra = collections.defaultdict(list)
        for chr_ in os.listdir(input):
            chr_dir = os.path.join(input, chr_)
            if os.path.isdir(chr_dir):
                for chunk_id in os.listdir(chr_dir):
                    chunk_dir = os.path.join(chr_dir, chunk_id)
                    if os.path.isdir(chunk_dir):
                        files = [os.path.join(chunk_dir, f) for f in os.listdir(chunk_dir) if f.endswith('.json')]
                        with Pool(cores) as pool:
                            for _, spectrum in pool.imap(f, files):
                                spectra[chr_ + '_' + chunk_id].append(spectrum)
    else:  # observed spectra
        spectra = {}
        for chr_ in os.listdir(input):
            chr_dir = os.path.join(input, chr_)
            if os.path.isdir(chr_dir):
                files = [os.path.join(chr_dir, f) for f in os.listdir(chr_dir) if f.endswith('.json')]
                with Pool(cores) as pool:
                    for chunk_id, spectrum in pool.imap(f, files):
                        spectra[chr_ + '_' + chunk_id] = spectrum

    np.savez(output, **{k: np.array(v) for k, v in spectra.items()})

    if table is not None:
        cols = ['chunk_id', 'maxp', 'power_maxp', 'snr_maxp']
        for i in range(low_p, high_p):
            cols += ['power_{0}p'.format(str(i)), 'snr_{0}p'.format(str(i))]
        df_dict = {k: [] for k in cols}
        x = np.linspace(low_p, high_p, num=resolution)

        if randomized:
            cols += ['item']
            df_dict.update({'item': []})

            for chunk_id, y_array in spectra.items():
                for i, y in enumerate(y_array):
                    df_dict['chunk_id'].append(chunk_id)
                    df_dict['item'].append(str(i))
                    maxp = x[np.argmax(y)]
                    df_dict['maxp'].append(maxp)
                    df_dict['power_maxp'].append(np.max(y))
                    snr_maxp = spectral.signal_to_noise(x, y, center=maxp)[0]
                    df_dict['snr_maxp'].append(snr_maxp)

                    for j in range(low_p, high_p):
                        df_dict['power_{0}p'.format(str(j))].append(y[spectral.closest_index(x, j)])
                        df_dict['snr_{0}p'.format(str(j))].append(spectral.signal_to_noise(x, y, center=j)[0])

        else:

            for chunk_id, y in spectra.items():
                df_dict['chunk_id'].append(chunk_id)
                maxp = x[np.argmax(y)]
                df_dict['maxp'].append(maxp)
                df_dict['power_maxp'].append(np.max(y))
                snr_maxp = spectral.signal_to_noise(x, y, center=maxp)[0]
                df_dict['snr_maxp'].append(snr_maxp)

                for j in range(low_p, high_p):
                    df_dict['power_{0}p'.format(str(j))].append(y[spectral.closest_index(x, j)])
                    df_dict['snr_{0}p'.format(str(j))].append(spectral.signal_to_noise(x, y, center=j)[0])

        df = pd.DataFrame(df_dict, columns=cols)
        df.to_csv(table, sep='\t', index=False)


if __name__ == '__main__':
    cli()
