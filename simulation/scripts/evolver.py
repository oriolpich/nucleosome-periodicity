import gzip
import pickle

import click
import numpy as np
from nucperiod import spectral
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()


def nucleosomes_positioning(index, dyad=200):
    """
    This will put nucleosomes in a random manner.
    """

    all_dyads = []

    # loop until all the sequence is filled
    while dyad + 73 < len(index):
        all_dyads.append(dyad)

        # add random linker distance
        random_linker = np.random.randint(30, 61)

        # select new dyad
        dyad = dyad + 73 + random_linker + 72

    return all_dyads


def return_signal(nucleosome_pos, seq):
    
    spline = robjects.r["smooth.spline"]
    xvector = np.arange(117)

    d_nucleotide = {0: 1 ,#'A',
                    1: 0, #'C',
                    2: 1, #'T',
                    3: 0, #'G'
                   }

    seqd = np.vectorize(d_nucleotide.get)(seq)
    array_nuc = []

    # select only the nucleosome positions
    for pos in nucleosome_pos:
        array_nuc.append(seqd[pos-58:pos+59])

    # do the stack
    center_nucleosome = np.sum(array_nuc, axis=0) / len(array_nuc)
    val = spline(xvector, center_nucleosome, )

    smoothed = list(val[1])

    x_s, y_s, snr_all, peak_random = spectral.compute(smoothed, low_t=0, high_t=115,
                           low_p=5, high_p=20, norm=True)

    snr10, peak = spectral.signal_to_noise(x_s, y_s, center=10.3)

    return peak_random, snr10, snr_all


def synthetic_signal_probability(P, prob, bounds=(-73, 74)):
    """Create synthetic signal probability"""
    x = np.arange(bounds[0], bounds[1])
    return np.cos(np.pi + 2 * np.pi * x / P)*prob


def evolve(seq, nucleosome_pos, mut_prob, mut_ype_prob, iterations, muts, length, probability):
    """ Method that evolves a sequence following three types of probability
    
    Args:
        iterations (int): Total number of iterations
        seq (np.array<np.int8>): a start sequence where 0=A,1=C,2=T,3=G (default is a random sequence)
        length (int): the length of the start sequence (if no sequence is provided)
        increase (np.array<np.int>): An array of positions that are 10% more probable
        decrease (np.array<np.int>): An array of positions that are 10% less probable        
        
    Returns:
        np.array: The sequence after the evolution
        np.array: The index where nucleosome dyads will be found
    
    """
    snr10_list = []
    snrall_list = []

    max_peak_list = []

    np_mut_prob = np.array([mut_prob[b] for b in "ACTG"], dtype=np.float64)
    np_mut_type_prob = np.array([[np.float64(mut_ype_prob[r].get(a, 0.0)) for a in "ACTG"] for r in "ACTG"])

    def base_mutation(ref):
        return np.random.choice(4, p=np_mut_type_prob[ref])

    mut_seq = np.vectorize(base_mutation)

    # define probability signal (wave-like)
    probability_signal = synthetic_signal_probability(10.3, probability, bounds=[-58, 59])

    probability_signal = 1 + probability_signal # this way it is already weighted!!

    # normalize the signal so that the sum is equal to one
    probability_signal = probability_signal/sum(probability_signal)

    # Start evolving the sequence
    evolved_sequences = []
    seq_original = seq.copy()
    ratio_vector = np.repeat(np.mean(probability_signal), len(seq))

    # ratio vector according to a wave affecting where the nucleosomes are located
    for pos in nucleosome_pos:
        ratio_vector[pos - 58:pos + 59] = probability_signal
    # normalize the ratios
    ratio_vector /= ratio_vector.sum()

    for i in range(iterations):
        
        # Initialize the sequence probability vector
        seq_prob = np_mut_prob[seq]
        
        # Apply the increase and decrease ratios
        seq_prob *= ratio_vector
        
        # Normalize the probability vector
        seq_prob /= seq_prob.sum()
        
        # Get mutated random positions
        mutated_positions = np.random.choice(length, p=seq_prob, size=muts)
        
        # Mutate the sequence
        alternates = mut_seq(seq[mutated_positions])
        seq[mutated_positions] = alternates
        
        # keep track of all changes
        evolved_sequences.append((mutated_positions, alternates))

        if i % 100 == 0:
            peak_random, snr10, snrall = return_signal(nucleosome_pos, seq)

            snr10_list.append(snr10)
            snrall_list.append(snrall)
            max_peak_list.append(peak_random)
        
    peak_random, snr10, snrall = return_signal(nucleosome_pos, seq)
    snr10_list.append(snr10)
    snrall_list.append(snrall)
    max_peak_list.append(peak_random)

    return seq_original, seq, nucleosome_pos, snr10_list, snrall_list, max_peak_list


def simulate(i, iterations, muts, length):

    np.random.seed(i)

    # Random initial sequence
    actg = np.array([0, 1, 2, 3], dtype=np.int8)
    seq = np.random.choice(actg, p=[0.25, 0.25, 0.25, 0.25], size=length, replace=True)

    index = np.arange(length)
    nucleosome_pos = nucleosomes_positioning(index, dyad=np.random.randint(150, 250))

    # Case 1

    # derived from Lynch paper in Nature EcoEvo
    prob_mut = {b: np.float64(p) for b, p in zip("ACTG", [0.06, 0.44, 0.06, 0.44])}

    prob_mut_type = {'A': {'C': 0.25, 'T': 0.25, 'G': 0.5},
                     'C': {'A': 0.125, 'T': 0.75, 'G': 0.125},
                     'T': {'G': 0.25, 'A': 0.25, 'C': 0.5},
                     'G': {'A': 0.75, 'T': 0.125, 'C': 0.125}}

    result1 = evolve(seq.copy(), nucleosome_pos, prob_mut, prob_mut_type, iterations=iterations, muts=muts,
                     length=length, probability=0.04)

    # Case 2

    result2 = evolve(seq.copy(), nucleosome_pos, prob_mut, prob_mut_type, iterations=iterations, muts=muts,
                     length=length, probability=0)

    # Case 3

    prob_mut = {b: np.float64(p) for b, p in zip("ACTG", [0.25, 0.25, 0.25, 0.25])}

    prob_mut_type = {'A': {'C': 1 / 3, 'T': 1 / 3, 'G': 1 / 3},
                     'C': {'A': 1 / 3, 'T': 1 / 3, 'G': 1 / 3},
                     'T': {'C': 1 / 3, 'A': 1 / 3, 'G': 1 / 3},
                     'G': {'C': 1 / 3, 'T': 1 / 3, 'A': 1 / 3}}

    result3 = evolve(seq.copy(), nucleosome_pos, prob_mut, prob_mut_type, iterations=iterations, muts=muts,
                     length=length, probability=0.04)

    return result1, result2, result3


@click.command()
@click.argument('output_prefix', metavar='<OUTPUT PREFIX>')
@click.option('-s', '--seed', type=click.INT, default=1, help='Random seed')
def cli(output_prefix, seed):
    """Perform 3 different evolutions of 6000 iterations, with 100 mutations each
    over a sequence of 1M nucleotides."""
    iterations = 6000
    muts = 100
    length = 1000000
    r1, r2, r3 = simulate(seed, iterations, muts, length)

    pickle.dump(r1, gzip.open('{}1_{}.pckl.gz'.format(output_prefix, seed), 'wb'))
    pickle.dump(r2, gzip.open('{}2_{}.pckl.gz'.format(output_prefix, seed), 'wb'))
    pickle.dump(r3, gzip.open('{}3_{}.pckl.gz'.format(output_prefix, seed), 'wb'))


if __name__ == '__main__':
    cli()
