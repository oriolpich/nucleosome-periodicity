import gzip
from os import path

import click

from nucperiod.utils import slicing_window


def get_polarized(input_folder, output_file):
    
    maf_files = [path.join(input_folder, 'chr{}.maf.gz'.format(name)) for name in list(range(1, 23)) + ['X', 'Y']]

    with gzip.open(output_file, 'wt') as out_polarized:
        for file in maf_files:
            with gzip.open(file, 'rt') as infile:
                for line in infile:
                    if line.startswith('#'):
                        continue
    
                    # Find alignment blocks
                    if 'a score' in line:
                        score = float(line.split('=')[1])
                        chrom, position, sequence_human, sequence_chimp, sequence_gorilla = None, None, None, None, None
                        if score > 0:
                            # read the whole block
                            next_line = next(infile)
                            while '\n' != next_line:
                                line_spl = next_line.split()
    
                                if line_spl[0] == 's':  # this is where the species and sequences are located
    
                                    specie = line_spl[1].split('.')[0]
    
                                    # iget the position where it is starting
                                    if specie == 'hg38':
                                        strand = line_spl[4]
                                        assert strand == '+'  # in theory al hg38 should be in the positive strand
                                        chrom = line_spl[1].split('.')[1]
                                        position = int(line_spl[2])
                                        sequence_human = line_spl[6].upper()
                                    elif specie == 'panTro4':
                                        sequence_chimp = line_spl[6].upper()
                                    elif specie == 'gorGor3':
                                        sequence_gorilla = line_spl[6].upper()
    
                                next_line = next(infile)
    
                            # process the block
                            if sequence_human is not None and sequence_gorilla is not None and sequence_chimp is not None:

                                human_offsets = {}
                                count = 0
                                for ix, i in enumerate(sequence_human):
                                    if i != '-':
                                        count += 1
                                    human_offsets[ix] = count
    
                                # get the triplets
                                triplets_chimp = list(slicing_window(sequence_chimp, 3))
                                triplets_gorilla = list(slicing_window(sequence_gorilla, 3))
                                triplets_human = list(slicing_window(sequence_human, 3))
    
                                for i, triplets in enumerate(zip(triplets_human, triplets_chimp, triplets_gorilla)):
                                    tri_human, tri_chimp, tri_gorilla = triplets
    
                                    if '-' in tri_human or '-' in tri_chimp or '-' in triplets_gorilla:
                                        continue
    
                                    pos = position + human_offsets[i + 1]
    
                                    ref = tri_chimp[1]
                                    alt = tri_human[1]
                                    if tri_chimp == tri_gorilla and \
                                            tri_chimp[0] == tri_human[0] and tri_chimp[2] == tri_human[2]:
                                        out_line = '{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos - 1, pos, ref, alt)
    
                                        if ref in 'CG':
                                            out_polarized.write(out_line)
                                        
                                        
@click.command()
@click.argument('input', metavar='<IN FOLDER>', type=click.Path(exists=True))
@click.argument('output', metavar='<OUT FILE>', type=click.Path())
def cli(input, output):
    """
    Parse the maf files to get polarized sites in human genome
    """
    get_polarized(input, output)


if __name__ == '__main__':
    cli()