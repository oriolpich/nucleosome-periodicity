
import gzip
import functools
import json
import itertools
import os
from collections import Counter
from multiprocessing import Pool

import bgdata
import click
import pandas as pd
from bgreference import refseq
from tqdm import tqdm


def read_file(f):
    """
    read the processed region file
    """
    df = pd.read_csv(f, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 's', 'e'])
    df['distance'] = df['e']-df['s']
    # return grouped
    return df.groupby(by='chr')


def kmers_generator(sequence, kmer):

    """
    Given the sequence, create tiplets
    :param sequence: a string
    :return: a list with the triplets
    """
    return list(slicing_window(sequence, kmer))


def slicing_window(seq, n=3):

    it = iter(seq)
    result = ''.join(itertools.islice(it, n))

    if len(result) == n:
        yield result

    for elem in it:
        result = result[1:] + elem
        yield result


def get_composition(group, kmer_len, species):
    """
    count the composition of each of the segments in a given k-mer and add it to the Counter
    """
    chrom, data = group
    total_triplets = Counter()

    for i, row in data.iterrows():
        try:
            seq = 0
            if kmer_len == 5:
                seq = refseq(species, chrom, int(row['s']) - 2, int(row['distance']) + 4)
            elif kmer_len == 7:
                seq = refseq(species, chrom, int(row['s']) - 3, int(row['distance']) + 6)
            elif kmer_len == 3:
                seq = refseq(species, chrom, int(row['s']) - 1, int(row['distance']) + 2)

            if len(seq)>0:
                total_triplets = total_triplets + Counter(kmers_generator(seq, kmer_len))
        except:
            continue

    return total_triplets


def get_full_composition(chrom, kmer_len, species):
    """
    count the composition of each of the segments and add it to the Counter
    """
    seq = refseq(species, chrom, 1, -1)
    return Counter(kmers_generator(seq, kmer_len))


@click.command()
@click.argument('species', metavar='<SPECIES>')
@click.argument('output', metavar='<OUT FILE>')
@click.option('-r', '--regions', type=click.Path(exists=True), default=None, help='Bed file with the regions of interest')
@click.option('--cores', type=click.INT, default=1, help='Number of cores to use for paralelization')
@click.option('--kmer', type=click.Choice(['3', '5', '7']), default='5', help='Size of the k-mer')
def cli(species, output, regions, cores, kmer):
    """
    Count the k-mers in a genome
    """
    kmer = int(kmer)

    count_composition = Counter()

    if regions is None:
        group = [n.replace('.txt', '') for n in os.listdir(bgdata.get_path('datasets', 'genomereference', species))
                 if not n.startswith('.') and not n.startswith('chrM') and not n.startswith('chr23') and not n.startswith('chr24')]
        
        f = functools.partial(get_full_composition, species=species, kmer_len=kmer)
    else:
        group = read_file(regions)
        f = functools.partial(get_composition, species=species, kmer_len=kmer)

    with Pool(int(cores)) as pool:
        for d in tqdm(pool.imap_unordered(f, group), total=len(group)):
            count_composition = count_composition + d

    with gzip.open(output, 'w') as fout:
        fout.write(json.dumps(dict(count_composition)).encode('utf-8'))


if __name__ == '__main__':
    cli()
