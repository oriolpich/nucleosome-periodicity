#=================================================================================
# Name          : expected_nucleosomes
# Author        : Oriol Pich
# Date          : 7/02/18
# Description   : It will calculat the number of expected mutations within the nucleosomes
#=================================================================================

import gzip
import json
import logging
import sys
from collections import defaultdict

import click
import numpy as np
from bgreference import refseq
import pandas as pd


# fix the random seed to always have the same results
np.random.seed(1234)


def compute_mutational_spectra(mutation_file, composition_file, kmers):
    """Count the distinct context (k-mer) for all mutations and divide by the genome counts
    of that k-mer"""
    df = pd.read_csv(mutation_file, sep='\t', names=['chr', 's', 'e', 'sample', 'kmer'])

    with gzip.open(composition_file, 'r') as fin:
        dic_composition = json.loads(fin.read().decode('utf-8'))

    if len(df):
        dic_mutations_type = df['kmer'].value_counts().to_dict()

        dic_probability_TN = None
        if kmers == 3:
            dic_probability_TN = defaultdict(lambda: defaultdict(dict))
        elif kmers == 5:
            dic_probability_TN = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))
        elif kmers == 7:
            dic_probability_TN = defaultdict(lambda: defaultdict(
                lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))))

        for k, count in dic_mutations_type.items():
            if 'N' not in k:
                if kmers == 3:
                    dic_probability_TN[k[0]][k[1]][k[2]] = count / dic_composition[k]
                elif kmers == 5:
                    dic_probability_TN[k[0]][k[1]][k[2]][k[3]][k[4]] = count / dic_composition[k]
                elif kmers == 7:
                    dic_probability_TN[k[0]][k[1]][k[2]][k[3]][k[4]][k[5]][k[6]] = count / dic_composition[k]
    else:
        logging.error('Not enough mutations to build spectra, exiting...')
        sys.exit(1)

    return dict(dic_probability_TN)


def kmers_generator_signature(sequence, spectra_dict, kmer):
    """
    Given the sequence, create the triplet probability. It feels odd but is the fastest way to compute
    :param sequence: given a sequence, create triplet and get the probability output
    """

    result = []
    seq_it = iter(sequence)

    if kmer == 3:
        a = next(seq_it)
        b = next(seq_it)
        for c in seq_it:
            try:
                signat_prob = spectra_dict[a][b][c]
            except:
                signat_prob = 0
            result.append(signat_prob)
            a, b = b, c
    elif kmer == 5:
        a = next(seq_it)
        b = next(seq_it)
        c = next(seq_it)
        d = next(seq_it)
        for e in seq_it:
            try:
                signat_prob = spectra_dict[a][b][c][d][e]
            except:
                signat_prob = 0
            result.append(signat_prob)
            a, b, c, d = b, c, d, e
    elif kmer == 7:
        a = next(seq_it)
        b = next(seq_it)
        c = next(seq_it)
        d = next(seq_it)
        e = next(seq_it)
        f = next(seq_it)
        for g in seq_it:
            try:
                signat_prob = spectra_dict[a][b][c][d][e][f][g]
            except:
                signat_prob = 0
            result.append(signat_prob)
            a, b, c, d, e, f = b, c, d, e, f, g

    return result


def get_sequence_window(nucid, species, window, kmers=3):
    """Get the sequence for a particular window (including flankig positions)"""
    chrom, start, end = nucid.split('_')

    # this is already one position before the real start
    start_real = int(start)-(window-1)/2
    t_win = 0

    if kmers == 3:
        t_win = window+2
        seq = refseq(species, chrom, int(start_real), t_win)
    elif kmers == 5:
        t_win = window+4
        seq = refseq(species, chrom, int(start_real) - 1, t_win)
    elif kmers == 7:
        t_win = window+6
        seq = refseq(species, chrom, int(start_real) - 2, t_win)

    # this is to avoid problems if we reach the end of the chromosome, unlikely in hg19 but can happen in yeast
    if len(seq) != t_win:
        diff = t_win-len(seq)
        seq = seq+('N'*diff)

    return seq


def bincount_app(rows, columns, n_rows, n_columns):

    # Get linear index equivalent
    lidx = (columns.max() + 1) * rows + columns

    # Use binned count on the linear indices
    return np.bincount(lidx, minlength=n_rows * n_columns).reshape(n_rows, n_columns)


def randomization_expected(possible_pos, mut_count, nrands, normalized_probability_vector,
                           final, len_win):

    # do all the randomizations requested at once
    matrix_rands = np.random.choice(possible_pos, size=(nrands, mut_count, ),
                                    p=list(normalized_probability_vector), replace=True)

    rows = np.repeat(range(nrands), mut_count)
    # create the column array
    cols = matrix_rands.reshape(1, mut_count * nrands)[0]

    final += bincount_app(rows, cols, nrands, len_win)

    return final


def get_expected_count_nucleosomes(species, mutation_file, composition_file, mutations_per_nucleosome_file, output_file, len_win, nrands, kmers=3):
    
    # create the mutational spectra dictionary to perform the randomizations
    spectra_dict = compute_mutational_spectra(mutation_file, composition_file, kmers)

    # create the matrix with all the randomizations to do
    matrix_rand = np.zeros((nrands, len_win))
    result = np.zeros(len_win)

    # this is for the randomization
    possible_pos = [i for i in range(len_win)]

    with gzip.open(mutations_per_nucleosome_file, 'rt') as infile:

        for i, line in enumerate(infile):

            nucid, mut_count = line.rstrip().split('\t')
            mut_count = int(mut_count)

            # get the sequence in the nucleosome fragment
            sequence = get_sequence_window(nucid, species, len_win, kmers)

            # get the vector with the probabilities after the mutational spectra correction
            vector_probabilities = kmers_generator_signature(sequence, spectra_dict, kmers)

            # Normalize the probabilities so that everything sum 1
            normalized_probability_vector = vector_probabilities / np.sum(vector_probabilities)

            # check if the vector is valid or not (probability mess?)
            if np.isnan(np.sum(normalized_probability_vector)):
                normalized_probability_vector = np.array([1/len_win]*len_win)

            # NOW WE WILL nrand times the location of the mutations nrands times
            if nrands > 0:
                matrix_rand = randomization_expected(possible_pos, mut_count, nrands,
                                                          normalized_probability_vector,
                                                     matrix_rand, len_win)

            # multiply all the probability by the observed count in the window
            expected_vector = normalized_probability_vector * int(mut_count)

            # here it will concatenate the vector
            result += expected_vector

    normalized_result = np.array(result)

    if nrands > 0:
        np.savez(output_file, freq=normalized_result, rand=matrix_rand)
    else:
        np.savez(output_file, freq=normalized_result)


@click.command()
@click.argument('species', metavar='<SPECIES>')
@click.argument('mutation_file', metavar='<MUTATIONS>')
@click.argument('composition_file', metavar='<COMPOSITION>')
@click.argument('mutations_per_nucleosome_file', metavar='<MUTATIONS PER NUCLEOSOME>')
@click.argument('output', metavar='<OUT FILE>')
@click.option('--size', 'len_win', type=click.INT, help='Size of the half window')
@click.option('--rands', 'nrands', type=click.INT, default=1000, help='Number of randomizations')
@click.option('--kmer', type=click.Choice(['3', '5', '7']), default='5', help='Size of the k-mer')
def cli(species, mutation_file, composition_file, mutations_per_nucleosome_file, output, len_win, nrands, kmer):
    """
    Count the k-mers in a genome
    """
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO, datefmt='%H:%M:%S')
    get_expected_count_nucleosomes(species, mutation_file, composition_file, mutations_per_nucleosome_file, output, len_win*2+1, nrands, int(kmer))


if __name__ == '__main__':
    cli()
