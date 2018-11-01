
import functools
import gzip
import itertools
import os
from multiprocessing import Pool

import click
import numpy as np


def get_from_fasta(fasta_file):
    """
    Args:
        fasta_file: file in .fa or .fa.gz format
    Returns:
        str: DNA sequence
    """
    with open(fasta_file, 'rt') as handle:
        fasta_str = ''
        for line in handle:
            if not line.startswith('>'):
                fasta_str += line.rstrip()
    return fasta_str


def markov_transition_matrix(template, k):

    dinucleotides = list(map(''.join, itertools.product(list('ACGT'), repeat=k)))
    matrix = np.zeros((4**k, 4**k))
    prev = template[:2]
    for i, a in enumerate(template[1:]):
        b = prev[1] + a
        if prev in dinucleotides and b in dinucleotides:
            matrix[dinucleotides.index(prev), dinucleotides.index(b)] += 1
        prev = b
    for j in range(4**k):
        s = sum(matrix[:, j])
        for i in range(4**k):
            matrix[i, j] /= s
    return matrix


def markov_stationary_state(transition_matrix, eps=1e-4, max_iter=1000):
    error = 1
    count = 0
    current = np.eye(transition_matrix.shape[0])[0, :]
    while (error > eps) and (count < max_iter):
        new = transition_matrix.dot(current)
        error = np.linalg.norm(current - new)
        current = new
    return current, error


def draw_from_markov_model(seed, trans_matrix, size):
    np.random.seed(seed=seed)
    n = trans_matrix.shape[0]
    dinucleotides = list(map(lambda x: x[0] + x[1], itertools.product(list('ACGT'), repeat=2)))
    stationary_state, error = markov_stationary_state(trans_matrix)
    initial = np.random.choice(n, p=list(stationary_state))
    draw = dinucleotides[initial]
    # prepare dict of random draws
    choices_dict = {d: list(np.random.choice(n,
                                  int(round((stationary_state[i] + 0.1) * size)),
                                  p=list(trans_matrix[:, i])
                                 )) for i, d in enumerate(dinucleotides)}
    new = draw
    for _ in range(size):
        b = choices_dict[new].pop()
        new = dinucleotides[b]
        draw += new[1]
    return draw


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--chunk', 'chunk_file', type=click.Path(), help='Chunk file')
@click.option('--output_folder', type=click.Path(), help='Folder where to save output files')
@click.option('--repeats', type=click.INT, default=100, help='Number of randomizations to perform')
@click.option('--cores', type=click.INT, default=1, help='Number of cores to use for paralelization')
def cli(chunk_file, output_folder, repeats, cores):
    """
    Generate a 16-state Markov
    """
    chunk = get_from_fasta(chunk_file)
    trans_matrix = markov_transition_matrix(chunk, 2)
    draw_func_partial = functools.partial(draw_from_markov_model, trans_matrix=trans_matrix, size=len(chunk))
    fmt = '{{:0{}d}}.txt.gz'.format(len(str(repeats-1)))
    with Pool(cores) as pool:
        for index, draw in enumerate(pool.imap(draw_func_partial, range(repeats))):
            out_file_name = fmt.format(index)
            with gzip.open(os.path.join(output_folder, out_file_name), 'wt') as g:
                g.write(draw)


if __name__ == '__main__':
    cli()
