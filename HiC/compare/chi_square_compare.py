import os
import re
from collections import namedtuple
import subprocess
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

import click
from cooler.api import Cooler
from scipy.stats import chisquare
import numpy as np

GenomeRange_ = namedtuple("GenomeRange", ["chr", "start", "end"])
class GenomeRange(GenomeRange_):
    def __str__(self):
        if (self.start is not None) and (self.end is not None):
            return "{}:{}-{}".format(self.chr, self.start, self.end)
        else:
            return self.chr
    
    def change_chromname(self):
        if self.chr.startswith("chr"):
            chr_ = self.chr.replace("chr", "")
            return GenomeRange(chr_, self.start, self.end)
        else:
            return GenomeRange("chr"+self.chr, self.start, self.end)

def region_str2genome_range(region):
    if '-' in region:
        chr_, s, e = re.split("[:-]", region)[:3]
        grange = GenomeRange(chr_, s, e)
    else:
        chr_ = region
        grange = GenomeRange(chr_, None, None)
    return grange



class MatrixSelector(object):
    """
    Selector for fetch the matrix from cool file.

    Parameters
    ----------
    cool : `cooler.api.Cooler`
        cool object.
    balance : bool
        balance matrix or not.
    process_func : callable, optional
        function for process the fetched array.
    """
    def __init__(self, cool, balance=True, process_func=None):
        self.cool = cool
        self.balance = balance
        self.process_func = process_func

    @property
    def chromsizes(self):
        return self.cool.chromsizes.to_dict()

    @property
    def binsize(self):
        return self.cool.binsize

    def binid_region2genome_range(self, binid_region):
        chr_, binid1, binid2 = binid_region
        assert chr_ in self.cool.chromnames
        resolution = self.binsize
        start = binid1 * resolution
        end = binid2 * resolution
        grange = GenomeRange(chr_, start, end)
        return grange

    def confirm_chromosome_name(self, region):
        """
        confirm region's chromosome in cool

        Parameters
        ----------
        region : str
        """
        grange = region_str2genome_range(region)
        chromnames = self.cool.chromnames
        if grange.chr not in chromnames:
            grange = grange.change_chromname()
            if grange.chr not in chromnames:
                raise IOError("chromosome {} not in cool file".format(grange.chr))
        return str(grange)

    def fetch(self, region1, region2=None):
        """
        Parameters
        ----------
        region1 : str
            like 'chr1:1000-2000'
        region2 : str
        """
        m = self.cool.matrix(balance=self.balance)
        region1 = self.confirm_chromosome_name(region1)
        if region2 is not None:
            region2 = self.confirm_chromosome_name(region2)
        arr =  m.fetch(region1, region2)
        if self.process_func is not None:
            arr = self.process_func(arr)
        return arr

    def fetch_by_bin(self, bin_region_1, bin_region_2):
        def convert_region(bin_region):
            if isinstance(bin_region, tuple):
                grange = self.binid_region2genome_range(bin_region)
                region = str(grange)
            else:
                chr_ = bin_region
                assert chr_ in self.cool.chromnames
                region = bin_region
            return region
        region1 = convert_region(bin_region_1)
        region2 = convert_region(bin_region_2)
        return self.fetch(region1, region2)


def read_bed(bed_path):
    """
    Read BED file.
    yield a series of genome ranges.
    """
    with open(bed_path) as f:
        for line in f:
            items = line.strip().split()
            chr_ = items[0]
            start = int(items[1])
            end = int(items[2])
            yield GenomeRange(chr_, start, end)


def chunking(block_iter, chunk_size):
    """
    group all blocks to several chunks,
    for parallel computing, 
    one process take one chunk as input.

    Parameters
    ----------
    block_iter : iterable
        A iterator yield many 'range blocks'.
    chunk_size : int
        How many blocks in one chunk.
    """
    chunk = []
    for idx, block in enumerate(block_iter):
        chunk.append(block)
        if (idx + 1) % chunk_size == 0:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

MIN_VAL = 1e-5

def process_region_chunk(range_chunk, mat_sel1, mat_sel2, gap_ratio_threshold):
    out_chunk = []
    for grange in range_chunk:
        mat1 = mat_sel1.fetch(str(grange))
        mat2 = mat_sel2.fetch(str(grange))
        nan1 = mat1[np.isnan(mat1)]
        nan2 = mat2[np.isnan(mat2)]
        gp_ratio_1 = nan1.shape[0] / (mat1.shape[0] * mat1.shape[1])
        gp_ratio_2 = nan2.shape[0] / (mat2.shape[0] * mat2.shape[1])
        if gp_ratio_1 > gap_ratio_threshold or gp_ratio_2 > gap_ratio_threshold:
            chi, p = np.nan, np.nan
        else:
            mat1 = mat1 / np.nansum(mat1)
            mat2 = mat2 / np.nansum(mat2)
            mat1[mat1 == 0] = MIN_VAL
            mat2[mat2 == 0] = MIN_VAL
            mat1[np.isnan(mat1)] = MIN_VAL
            mat2[np.isnan(mat2)] = MIN_VAL
            mat1 = mat1.reshape(-1)
            mat2 = mat2.reshape(-1)
            chi, p = chisquare(mat2, f_exp=mat1)
        res = (grange.chr, grange.start, grange.end, chi, p, gp_ratio_1, gp_ratio_2)
        out_chunk.append(res)
    return out_chunk


def write_result(out_chunk, output, mode="a"):
    with open(output, mode=mode) as f:
        for res in out_chunk:
            outline = "\t".join([str(i) for i in res]) + "\n"
            f.write(outline)


@click.command()
@click.argument("cool-uri-1")
@click.argument("cool-uri-2")
@click.argument("bed-path")
@click.argument("output")
@click.option("--inner-window", "-i",
    type=int, default=3,
    show_default=True,
    help="The size of inner window, "
         "interactions within inner window will be removed "
         "as self-ligation interaction. "
         "The unit is number of bins.")
@click.option("--balance/--no-balance",
    default=True,
    show_default=True,
    help="Use balanced matrix or not.")
@click.option("--gap-ratio", "-g",
    default=0.2,
    show_default=True,
    help="The gap ratio threshold, "
         "only gap ratio less equal than this will be keeped. ")
@click.option("--processes", "-p",
    type=int, default=1,
    show_default=True,
    help="Number of processes.")
@click.option("--chunk-size", "-s",
    type=int, default=5000,
    show_default=True,
    help="How many blocks in one chunk, for parallel processing.")
def chi_square(cool_uri_1, cool_uri_2, bed_path, output, inner_window, balance, gap_ratio, processes, chunk_size):
    """
    \b
    Args
    ----
    cool_uri_1 : str
        URI of input cool container 1.
    cool_uri_2 : str
        URI of input cool container 2.
    bed_path : str
        Path to input BED.
    output : str
        Path to output BEDGRAPH file.
    """
    c1 = Cooler(cool_uri_1)
    matrix_selector_1 = MatrixSelector(c1, balance=balance)
    c2 = Cooler(cool_uri_2)
    matrix_selector_2 = MatrixSelector(c2, balance=balance)
    
    regions = read_bed(bed_path)
    chunks = chunking(regions, chunk_size)
    if os.path.exists(output):
        subprocess.check_call(['rm', output])
    with ProcessPoolExecutor(max_workers=processes) as excuter:
        map_ = map if processes == 1 else excuter.map
        for out_chunk in map_(process_region_chunk, chunks, repeat(matrix_selector_1), repeat(matrix_selector_2), repeat(gap_ratio)):
            write_result(out_chunk, output, mode='a')


if __name__ == "__main__":
    eval("chi_square()")