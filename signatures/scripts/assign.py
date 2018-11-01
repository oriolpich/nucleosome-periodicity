
import gzip
import logging
import os

import click
import pandas as pd
from collections import defaultdict


def get_real_mut(df):

    reverse = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    spl = df['Mutation_type'].split('[')
    first = spl[0]
    spl = df['Mutation_type'].split('[')[1].split('>')[0]
    second = spl
    spl = df['Mutation_type'].split('>')[1].split(']')[1]
    third = spl

    spl = df['Mutation_type'].split(']')[0].split('>')[1]
    fourth = spl

    mut = '{}{}{}_{}'.format(first, second, third, fourth)
    reversed = '{}{}{}_{}'.format(reverse[third], reverse[second], reverse[first], reverse[fourth])

    df['m'] = mut
    df['rev'] = reversed

    return df


def load_cosmic(cosmic_file):
    W = pd.read_csv(cosmic_file, sep='\t')
    W.dropna(axis=1, inplace=True)
    renames = {c: '.'.join(c.split()) for c in W.columns if c.startswith('Signature')}
    W.set_index('Somatic Mutation Type', inplace=True)
    W.rename(columns=renames, inplace=True)
    W.drop(columns=['Substitution Type', 'Trinucleotide'], inplace=True)

    return W


def assign_probability(inputfile, weights, cosmic_file, output_dir):

    W = load_cosmic(cosmic_file)

    H = pd.read_csv(weights, header=0, sep="\t")

    if 'signatures_weight' in weights:
        SAMPLE = 'sample_id'
        CHR_PREFIX = 'chr'
        ALT_INCLUDED = True
    elif 'exposures' in weights:
        SAMPLE = 'Sample'
        CHR_PREFIX = ''
        ALT_INCLUDED = False
        df_original = pd.read_csv(inputfile, sep='\t', names=['sample', 'chr', 'pos', 'ref', 'alt', 'triplet'])
        number_muts = df_original['sample'].value_counts().to_dict()
        H['mutation_count'] = H[SAMPLE].map(number_muts)
        # fix the columns
        Hcols = H.columns
        new_cols = [col.replace(' ', '.') for col in Hcols]
        H.columns = new_cols
    else:
        raise ValueError('Unknown project')

    # go over each sample in H matrix and compute the probability for each mutation type in a tri-nuc context

    frames = []  # to collect results sample wise
    flag = 0
    for idx, row in H.iterrows():  # go over each sample
        sample = row[SAMPLE]
        sig_dic = {}
        allsigs = []
        # get the exposure (i.e total number of mutations belong to each signature) value for the particular sample from H matrix
        for col in H.columns:
            if col not in [SAMPLE, 'SSE', 'mutation_count']:
                sig_dic[col] = row[col] * row['mutation_count']  # save the exposure value in a dictionary per signature name
                allsigs.append(col)  # save the signature names

        # multiple the exposure (from H) with the W matrix
        a = W.copy()  # take a copy of the W matrix (which is the extracted signatures - not sample specific)
        for sig in allsigs:
            a[sig] *= sig_dic[sig]  # multiply the signature columns with the corresponding signature exposure in that particular sample

        # compute the row sum for normalization (i.e sum of values across signature for each mutation/context type)
        a['row_sum'] = a[allsigs].sum(axis=1)

        # normalize the row values with the row sum to driver
        # the probabilities for different signatures for each mutation type
        new = a[allsigs].div(a['row_sum'], axis=0)[allsigs]

        # add info columns
        new['Mutation_type'] = new.index
        new['Sample'] = sample

        # sort the columns
        columns = ['Sample', 'Mutation_type'] + allsigs

        new = new[columns]

        # save the results for each samples in a dataframe
        if flag == 0:
            frames = [new]
            flag += 1
        else:
            frames.append(new)

    results_new = pd.concat(frames)
    results_new.dropna(inplace=True)

    selected = [col for col in results_new.columns if col.startswith('Signature')]
    res = results_new[selected].idxmax(axis=1)

    df = results_new[['Sample', 'Mutation_type']].copy()
    df['Max'] = res

    # get mutation labels
    df = df.apply(get_real_mut, axis=1)

    # dictionary where we will store where everything goes
    dic_sig = defaultdict(lambda: defaultdict(str))
    all_sigs_found = set()

    # get the maximum of each row
    for i, row in df.iterrows():
        dic_sig[row['Sample']][row['m']] = row['Max']
        dic_sig[row['Sample']][row['rev']] = row['Max']
        all_sigs_found.add(row['Max'])

    # assign each mutation to each signature

    # open all files
    all_files_out = {s: gzip.open(os.path.join(output_dir, '{}.tsv.gz'.format(s.replace('.', '_'))), 'wt')
                     for s in all_sigs_found}

    missing = defaultdict(int)
    # assign each mutation
    with gzip.open(inputfile, 'rt') as infile:

        for line in infile:

            line_spl = line.rstrip().split('\t')
            donor = line_spl[0]
            ref = line_spl[3]
            alt = line_spl[4]

            if ALT_INCLUDED:
                triplet = line_spl[-1]
            else:
                triplet = '{}_{}'.format(line_spl[-1], alt)

            if ('N' not in triplet) and (line_spl[1] != 'M'):

                out = '{}{}\t{}\t{}\t{}\t{}\n'.format(CHR_PREFIX, line_spl[1], line_spl[2], ref, alt, donor)

                try:
                    all_files_out[dic_sig[donor][triplet]].write(out)
                except KeyError:
                    missing[donor] += 1
                    logging.debug('Unable to write %s' % out)

    for file in all_files_out.values():
        file.close()

    for donor, counts in missing.items():
        logging.warning('Unable to write %s muts for sample %s' % (counts, donor))


@click.command()
@click.argument('input', metavar='<IN FILE>', type=click.Path(exists=True))
@click.argument('weights', metavar='<WEIGHTS FILE>', type=click.Path(exists=True))
@click.argument('signatures', metavar='<SIGNATURES FILE>', type=click.Path(exists=True))
@click.argument('output', metavar='<OUT FOLDER>', type=click.Path(exists=True))
@click.option('--verbose', '-v', is_flag=True, default=False, help='Increase verbosity')
def cli(input, weights, signatures, output, verbose):
    """
    Assign each mutation to the most probable signature it comes from.

    The assigment of signatures is done by either DeconstructSigs
    or SigFit R packages
    """
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                        level=logging.DEBUG if verbose else logging.INFO, datefmt='%H:%M:%S')
    assign_probability(input, weights, signatures, output)


if __name__ == '__main__':
    cli()
