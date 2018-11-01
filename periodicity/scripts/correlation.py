
import glob
import gzip
import json
import multiprocessing
import os

import click


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


def is_head(chunk, motif):
    """
    Args:
        chunk: str
        motif: set of str
    Returns:
        bool: whether some motif in motif_set is head of chunk
    """
    for s in motif:
        if chunk.startswith(s):
            return True
    return False


def pair_count(chunk, motif, dist_bound=120):
    """
    Args:
        chunk: str: DNA string
        motif: set of str with same length = 2
        dist_bound: int: upper bound for origin-end distance
    Returns:
        dict with:
            key: 'pair_count'; value: array: entries at index d represent pair counts at distance d
            key: 'motif_position'; value: array: indices in chunk where the motif is found:
                                   first position is index, subsequent positions are increments
            key: 'motif_count'; value: int: motif count
            key: 'chunk_len'; value: int: length of chunk
    """
    motif_pos = []
    for i in range(len(chunk)):
        if is_head(chunk[i: i + 2], motif):
            motif_pos.append(i)
    pc = [0 for _ in range(dist_bound)]
    for i, _ in enumerate(motif_pos):
        j = i + 1
        if j >= len(motif_pos):
            continue
        a = motif_pos[j] - motif_pos[i]
        while a < dist_bound:
            pc[a] += 1
            j += 1
            if j >= len(motif_pos):
                break
            a = motif_pos[j] - motif_pos[i]
    return {'pair_count': pc, 'motif_count': len(motif_pos), 'chunk_len': len(chunk)}


def main(input_file, output_file):
    motif = {'AA', 'TA', 'AT', 'TT'}
    if input_file.endswith('.fa'):
        chunk = get_from_fasta(input_file)
    else:  # assume gzip file
        with gzip.open(input_file, 'rt') as f:
            chunk = f.read()
    res = pair_count(chunk, motif)
    with open(output_file, 'w+') as f:
        json.dump(res, f)


def mapper(t):
    return main(*t)


@click.command()
@click.argument('input', metavar='<DIR with FASTA>', type=click.Path(exists=True))
@click.argument('output', metavar='<OUTPUT FOLDER>', type=click.Path())
@click.option('--cores', type=click.INT, default=1, help='Number of cores to use for paralelization')
@click.option('--wildcard', default=None, help='Wildcard to pass to the glob module')
def cli(input, output, cores, wildcard):
    """Compute motifs"""
    if wildcard is None:
        in_files = [os.path.join(input, f) for f in os.listdir(input)]
    else:
        in_files = [f for f in glob.glob(os.path.join(input, wildcard))]

    names = []
    for file in in_files:
        name = os.path.basename(file)
        if file.endswith('.fa'):
            names.append(name.rsplit('.', 1)[0])
        else:  # assume gzip
            names.append(name.rsplit('.', 2)[0])

    out_files = [os.path.join(output, name+'.json') for name in names]

    with multiprocessing.Pool(cores) as pool:
        for _ in pool.imap(mapper, zip(in_files, out_files)):
            pass


if __name__ == '__main__':
    cli()
