import glob
from os import path

import click
import numpy as np
import pandas as pd


def join_signatures(infolder, outfolder):

    # get all the 30 cosmic signatures
    cosmic_signatures = np.arange(1, 31)

    # create the name of the cosmic signatures
    name_cosmic_signatures = ['Signature_{}.tsv.gz'.format(sig) for sig in cosmic_signatures]

    # iterate over each signature and concatenate all the dfs
    for signature in name_cosmic_signatures:
        path_ = path.join(infolder, '*', signature)
        files = glob.glob(path_)

        if len(files) > 0:
            all_dfs = []
            for f in files:
                df = pd.read_csv(f, sep ='\t', names=['chr', 'pos', 'ref','alt', 'sample'])
                all_dfs.append(df)

            sig_df = pd.concat(all_dfs)

            # save the output to the specified directory
            outf = path.join(outfolder, signature)
            sig_df.to_csv(outf, sep='\t', index=False, header=False, compression='gzip')


@click.command()
@click.argument('input', metavar='<IN FOLDER>', type=click.Path(exists=True))
@click.argument('output', metavar='<OUT FOLDER>', type=click.Path(exists=True))
def cli(input, output):
    """
    Concatenate all signatures files that correspond to the same signatures.

    The files are search in inputfolder/*/Signature_X.tsv.gz
    """
    join_signatures(input, output)


if __name__ == '__main__':
    cli()
