import glob
import gzip
import os

import click
import pandas as pd


def split_file_samples(file, output_folder):
    keep_track = []

    df = pd.read_csv(file, sep='\t', names=['chr', 'pos', 'ref', 'alt', 'sample'])
    for sample, data in df.groupby(by='sample'):
        outfile = '{}/{}.tsv.gz'.format(output_folder, sample)
        data.to_csv(outfile, sep='\t', index=False, header=False, compression='gzip')
        out = '{}\t{}\t{}\n'.format(sample, os.path.basename(file), len(data))
        keep_track.append(out)
    return keep_track


def split_samples(project_folder, output_folder):
    keep_track = []

    for file in glob.iglob('{}/*.tsv.gz'.format(project_folder)):
        basename = os.path.basename(file)
        ctype_folder = os.path.join(output_folder, basename.replace('.tsv.gz', ''))
        os.makedirs(ctype_folder, exist_ok=True)
        keep_track += split_file_samples(file, ctype_folder)

    # write the info in a file
    with gzip.open('{}/sample_tracksheet.tsv.gz'.format(output_folder), 'wt') as outfile:
        for line in keep_track:
            outfile.write(line)


@click.command()
@click.argument('input_folder', metavar='<IN FOLDER>', type=click.Path(exists=True))
@click.argument('output_folder', metavar='<OUT FOLDER>', type=click.Path())
def cli(input_folder, output_folder):
    """
    Find all .tsv.gz file in a folder and split
    split those files into one file per sample.

    Generates a summary as ``sample_tracksheet.tsv.gz``
    """
    os.makedirs(output_folder, exist_ok=True)
    split_samples(input_folder, output_folder)


if __name__ == '__main__':
    cli()
