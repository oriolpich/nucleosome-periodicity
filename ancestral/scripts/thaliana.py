import gzip
from os import path

import click


def parse_mafs(folder, output, species_of_interest):
    name = path.basename(folder)
    maf_files = [path.join(folder, '{}.chr{}.maf'.format(name, chr_)) for chr_ in range(1, 6)]

    with gzip.open(output, 'wt') as outf:
        for file in maf_files:
            with open(file) as inf:
                for line in inf:
                    # remove beginning of the file
                    if line.startswith('#'):
                        continue

                    if 'a#' in line:
                        chrom, position, sequence_thaliana, sequence_interest = None, None, None, None
                        next_line = next(inf)
                        # read block
                        while '\n' != next_line:  # every block ends with an empty line
                            line_spl = next_line.split()
                            if line_spl[0] == 's':  # this is where the species and sequences are located

                                species = line_spl[1].split('.')[0]

                                # if desired specie, get the position where it is starting
                                if species == 'arabidopsis_thaliana':
                                    strand = line_spl[4]
                                    assert strand == '+'
                                    chrom = line_spl[1].split('.')[1]
                                    position = int(line_spl[2])
                                    sequence_thaliana = line_spl[6].upper()
                                elif species == species_of_interest:
                                    sequence_interest = line_spl[6].upper()

                            next_line = next(inf)

                        if sequence_thaliana is not None and sequence_interest is not None:
                            # process the block

                            count = 0

                            # get all real_positions and changes to compare
                            for ix, nuc in enumerate(sequence_thaliana):
                                if nuc != '-':
                                    count += 1
                                if ix > 0 & ix < len(sequence_thaliana) - 1:
                                    tri_interest = ''.join(sequence_interest[ix - 1:ix + 2])
                                    tri_thaliana = ''.join(sequence_thaliana[ix - 1:ix + 2])
                                    out = '{}\t{}\t{}\t{}\t{}\n'.format(chrom, position + count - 1, position + count,
                                                                        tri_interest, tri_thaliana)
                                    if '-' not in out:
                                        outf.write(out)


def get_polarized(in_file, out_file):
    with gzip.open(in_file, 'rt') as infile, \
            gzip.open(out_file, 'wt') as outfile:
        for line in infile:
            line_spl = line.rstrip().split('\t')
            nucls_outgroup = line_spl[4]
            thaliana_nuc = line_spl[5]

            ref = nucls_outgroup[1]
            alt = thaliana_nuc[1]

            if len(nucls_outgroup) == 3 and len(thaliana_nuc) == 3 and nucls_outgroup[0] == thaliana_nuc[0] and \
                    nucls_outgroup[-1] == thaliana_nuc[-1] and ref in 'CG':
                out = 'chr{}\t{}\t{}\t{}\t{}\n'.format(line_spl[0], line_spl[1], line_spl[2], ref, alt)
                outfile.write(out)


@click.group()
def cli():
    """
    Utilities for processing thaliana data
    """
    pass


@cli.command()
@click.argument('input', metavar='<IN FOLDER>', type=click.Path(exists=True))
@click.argument('species', metavar='<SPECIES>')
@click.argument('output', metavar='<OUT FILE>', type=click.Path())
def parse(input, species, output):
    """
    Parse maf files associated with a particular species
    """
    parse_mafs(input, output, species)


@cli.command()
@click.argument('input', metavar='<IN FOLDER>', type=click.Path(exists=True))
@click.argument('output', metavar='<OUT FILE>', type=click.Path())
def find(input, output):
    """
    Get the polarized sites
    """
    get_polarized(input, output)


if __name__ == '__main__':
    cli()
