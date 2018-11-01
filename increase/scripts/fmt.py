
import click
import pandas as pd
from bgreference import refseq


header = ['CHR', 'POS', 'REF', 'ALT', 'SAMPLE']


def get_mutation(row, genome, size):
    ref = row['REF']
    chr = row['CHR']
    pos = int(row['POS'])
    try:
        kmer = refseq(genome, chr, pos - 2, 5)
    except (ValueError, RuntimeError):  # out of chrom size, unknown chr
        return None
    if len(kmer) == 5 and (ref is '-' or kmer[2] == ref) and all(n in 'ACGT' for n in kmer):
        if size == 5:
            return kmer
        else:  # 3-mer by default
            return kmer[1:4]
    else:
        return None


def main(input, output, genome, kmer):
    df = pd.read_csv(input, sep='\t', header=None)
    header_size = len(df.columns)
    df.rename(columns={i: n for i,n in enumerate(header[:header_size])}, inplace=True)
    df['MUT'] = df.apply(get_mutation, axis=1, genome=genome, size=kmer)
    df = df[df['CHR'] != 'M']
    df = df[df['CHR'] != 'chrM']
    df = df[~df['MUT'].isnull()]
    df['POS-1'] = df['POS'] - 1
    df.to_csv(output, header=False, compression='gzip', sep='\t', index=False,
              columns=['CHR', 'POS-1', 'POS', 'SAMPLE', 'MUT'])


@click.command()
@click.argument('input', metavar='<MUT FILE>')
@click.argument('output', metavar='<OUT FILE>')
@click.argument('genome', metavar='<GENOME>')
@click.argument('size', metavar='<K-MER size>', type=click.INT)
def cli(input, output, genome, size):
    main(input, output, genome, size)


if __name__ == '__main__':
    cli()
