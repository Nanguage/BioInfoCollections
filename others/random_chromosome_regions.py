import random

import numpy as np
import click


@click.group()
def random_chromosome_regions():
    pass


class RandomChromGen(object):
    def __init__(self, chromsizes):
        self.chromsizes = chromsizes
        chroms, lengths = list(zip(*chromsizes.items()))
        self.chroms = list(chroms)
        self.lengths = lengths
        tol_len = sum(lengths)
        self.pk = [l/tol_len for l in lengths]
    
    def random_chrom(self):
        return np.random.choice(self.chroms, p=self.pk)


def load_chromsizes(path):
    chromsizes = {}
    with open(path) as f:
        for line in f:
            items = line.strip().split()
            chromsizes[items[0]] = int(items[1])
    return chromsizes


def random_fixed_length_region(chrom_length, region_length):
    start = np.random.randint(chrom_length - region_length)
    end = start + region_length
    return start, end


@random_chromosome_regions.command()
@click.argument("chromosome-sizes")
@click.argument("output")
@click.option("--length", "-l",
    default=1000,
    show_default=True,
    help="The length of regions.")
@click.option("--number", "-N",
    default=100,
    show_default=True,
    help="Number of random regions.")
@click.option("--prob-forward", "-p",
    default=0.5, type=float,
    show_default=True,
    help="Probability of forward strand.")
@click.option("--seed", "-s",
    default=1,
    show_default=True,
    help="Random seed.")
def fixed_length(chromosome_sizes, output, length, number, prob_forward, seed):
    """
    Random select a series of chromosome regions in fixed length,
    output in BED6 format.

    \b
    Args
    ----
    chromosome_sizes : str
        Path to chromosome size file.
        Chromosome size file is a 2 columns tab delimited file,
        first columns is chromosome name, second column is it's length.
    output : str
        Path to output file.
    """
    np.random.seed(seed=seed)
    chromsizes = load_chromsizes(chromosome_sizes)
    chrom_gen = RandomChromGen(chromsizes)
    with open(output, 'w') as f:
        for _ in range(number):
            chr_ = chrom_gen.random_chrom()
            start, end = random_fixed_length_region(chromsizes[chr_], length)
            strand = np.random.choice(['+', '-'], p=[prob_forward, 1-prob_forward])
            out = [chr_, start, end, '.', 0, strand]
            outline = "\t".join([str(i) for i in out])
            f.write(outline+"\n")


if __name__ == "__main__":
    eval("random_chromosome_regions()")
