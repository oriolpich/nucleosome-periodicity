import sys

import click
import pandas as pd
from bgreference import hg19


def get_mutation_deconstructsigs(row):
    ref = row['REF']
    alt = row['ALT']
    chr = row['CHR']
    pos = int(row['POS'])
    triplet = hg19(chr, pos - 1, 3)
    if triplet[1] == ref and 'N' not in triplet:
        return '_'.join([triplet, alt])
    else:
        return None


def fmt_deconstructsigs(infile, outfile, min_samples):
    df = pd.read_csv(infile, sep='\t', header=None, names=['CHR', 'POS', 'REF', 'ALT', 'SAMPLE'])
    df['CHR'] = df['CHR'].str.replace('chr', '')
    df['MUT'] = df.apply(get_mutation_deconstructsigs, axis=1)
    df = df[df['CHR'] != 'M']
    df = df[~df['MUT'].isnull()]
    df.sort_values(by=['SAMPLE', 'CHR'])
    # Filter samples with less than X mutations
    sample_counts = df['SAMPLE'].value_counts()
    df = df[df['SAMPLE'].isin(sample_counts.index[sample_counts > min_samples])]

    if df.empty:
        sys.exit(3)
    else:
        df.to_csv(outfile, header=False, compression='gzip', sep='\t',
                  index=False,
                  columns=['SAMPLE', 'CHR', 'POS', 'REF', 'ALT', 'MUT'])


def get_mutation_sigfit(row):
    ref = row['REF']
    chr = row['CHR']
    pos = int(row['POS'])
    triplet = hg19(chr, pos - 1, 3)
    if triplet[1] == ref and 'N' not in triplet:
        return triplet
    else:
        return None


def fmt_sigfit(infile, outfile, min_samples):
    df = pd.read_csv(infile, sep='\t', header=None, names=['CHR', 'POS', 'REF', 'ALT', 'SAMPLE'])
    df['MUT'] = df.apply(get_mutation_sigfit, axis=1)
    df = df[~df['MUT'].isnull()]
    df.sort_values(by=['SAMPLE', 'MUT', 'CHR'])
    # Filter samples with less than X mutations
    sample_counts = df['SAMPLE'].value_counts()
    df = df[df['SAMPLE'].isin(sample_counts.index[sample_counts > min_samples])]

    if df.empty:
        sys.exit(3)
    else:
        df.to_csv(outfile, header=False, compression='gzip', sep='\t', index=False,
                  columns=['SAMPLE', 'CHR', 'POS', 'REF', 'ALT', 'MUT'])


@click.command()
@click.argument('input', metavar='<IN FILE>', type=click.Path(exists=True))
@click.argument('output', metavar='<OUT FILE>', type=click.Path())
@click.option('--min', 'min_samples', default=50)
def cli(input, output, min_samples):
    """
    Format a mutations file to be able to be used by DeconstructSigs or Sigfit
    """
    if 'deconstructsigs' in output:
        fmt_deconstructsigs(input, output, min_samples)
    elif 'sigfit' in output:
        fmt_sigfit(input, output, min_samples)
    else:
        raise ValueError('Unknown project')


if __name__ == '__main__':
    cli()
