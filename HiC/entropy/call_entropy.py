"""
"""

import os
import re
from itertools import repeat, tee
import operator
from concurrent.futures import ProcessPoolExecutor
import subprocess
from collections import namedtuple

import click
from cooler.api import Cooler
import numpy as np
import scipy.stats
from sklearn.linear_model import LinearRegression


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


BedGraph = namedtuple("BedGraph", ['chrom', 'start', 'end', 'value'])


def tiling(bin_num, window_size, overlap_size):
    """
    yield a series of 'range block',
    a 'range block' means a start position and a end position,
    in bin ticks (one tick one bin).

                                  -----  block N
                    ...
                                               ...
          -----
       -----                                   2
    -----                                block 1
    |---------------------------------|
    0                            (bin_num - 1)

    Parameters
    ----------
    bin_num : int
        Total bin number
    window_size : int
    overlap : int
    """
    start = 0
    end = start + window_size
    step = window_size - overlap_size
    while 1:
        if end > (bin_num - 1):
            end = bin_num - 1
            yield (start, end)
            break
        elif end == (bin_num - 1):
            yield (start, end)
            break
        else:
            yield (start, end)
            start += step
            end = start + window_size


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


def chromosome_chunks(chromsizes, window_size, overlap, chunk_size):
    """
    yield chunks with chromosome name 

    Parameters
    ----------
    chromsizes : dict
        mapping from chromosome name to it's size, unit in bins.
    window_size : int
    overlap : int
    chunk_size : int
    """
    for chr_, bin_num in chromsizes.items():
        blocks = tiling(bin_num, window_size, overlap)
        chunks = chunking(blocks, chunk_size)
        for ck in chunks:
            yield (chr_, ck)


def eliminate_inner(mat, start, end, inner_window):
    """
    eliminate inner window of fetched matrix in loci entropy calling.
    """
    flank = (inner_window - 1) // 2
    inner_start = max(0, start-flank)
    inner_end = min(end+flank, mat.shape[1])
    inner_mat = mat[:, inner_start:inner_end]
    s = inner_mat.shape
    mask = sum([np.eye(*s, k=i) for i in np.arange(-flank, flank+1)])
    inner_mat[mask == 1] = np.nan
    return mat


def cut_outer(mat, start, end, outer_window):
    """
    Exclude the interactions, outside the outer window.

    Return
    ------
    res : `numpy.ndarray`
        Cutted matrix.
    c : int
        Relative matrix center index.
    """
    center = (start + end) // 2
    flank = (outer_window - 1) // 2
    if center - flank < 0 and center + flank > mat.shape[1]:
        res = mat
        c = res.shape[1] // 2 + 1
    elif center - flank < 0:
        res = mat[:, :outer_window]
        c = center
    elif center + flank > mat.shape[1]:
        res = mat[:, -outer_window:]
        c = center - (mat.shape[1] - outer_window)
    else:
        res = mat[:, center-flank:center+flank+1]
        c = res.shape[1] // 2 + 1
    return res, c


def subtract_arr_expect(arr, relative_center):
    """
    Estimate the expect value and subtract it.
    """
    MIN_NUM = 3
    up_part = arr[:relative_center]
    up_part = up_part[::-1]
    down_part = arr[relative_center:]
    if up_part.shape[0] >= MIN_NUM and down_part.shape[0] >= MIN_NUM:
        expect_up = estimate_expect_1d(up_part)
        expect_down = estimate_expect_1d(down_part)
        up_part_ = (up_part - expect_up)[1:]
        down_part_ = (down_part - expect_down)[1:]
        res = np.concatenate((up_part_[::-1], np.asarray([0]), down_part_))
    elif up_part.shape[0] < MIN_NUM and down_part.shape[0] >= MIN_NUM:
        expect_down = estimate_expect_1d(down_part)
        down_part_ = (down_part - expect_down)[1:]
        res = np.concatenate( (np.zeros(MIN_NUM), down_part_) )
    elif up_part.shape[0] >= MIN_NUM and down_part.shape[0] < MIN_NUM:
        expect_up = estimate_expect_1d(up_part)
        up_part_ = (up_part - expect_up)[1:]
        res = np.concatenate( (up_part_, np.zeros(MIN_NUM)) )
    else:
        raise ValueError("Input array too small")
    return res


def estimate_expect_1d(arr):
    xs = np.arange(arr.shape[0])
    xs = xs[arr != 0]
    xs = np.log2(xs + 1)
    ys = np.log2(arr[arr!=0])
    xs = xs[~np.isnan(ys)]
    ys = ys[~np.isnan(ys)]
    if xs.shape[0] > 1:
        lr = LinearRegression()
        lr.fit(X=xs.reshape(-1,1), y=ys.reshape(-1,1))
        a, b = lr.coef_[0][0], lr.intercept_[0]
        def f(d):
            return (2 ** b) * (d ** a)
        expect = f(np.arange(1, arr.shape[0]+1))
        return expect
    else:
        return np.zeros(arr.shape[0])


def relu(arr):
    return np.where(arr > 0, arr, 0)


def sliding_cal_entropy(selector, chrom, chunk, inner_window, outer_window, subtract_expect, non_nan_threshold=0.6):
    """
    iterate over the chunk, calculate a value.

    Parameter
    ---------
    selector : `MatrixSelector`
        selector for fetch the contact matrix.
    chrom : str
        chromosome name.
    chunk : iterable
        blocks
    inner_window : int
        inner window size.
    outer_window : int
        outer window size.
    minus_except : bool
        minus except value or not.
    non_nan_threshold : float
        non-NaN value threshold, if less than this will be droped.
    """
    start, end = chunk[0][0], chunk[-1][-1]
    matrix = selector.fetch_by_bin((chrom, start, end), chrom)
    if inner_window > 0:
        matrix = eliminate_inner(matrix, start, end, inner_window)
    for s, e in chunk:
        s_, e_ = s - start, e - start
        m = matrix[s_:e_, :]
        if outer_window > 0:
            m, c = cut_outer(m, s, e, outer_window)
        arr = np.nanmean(m, axis=0)
        non_nan_rate = arr[~np.isnan(arr)].shape[0] / arr.shape[0]
        if non_nan_rate < non_nan_threshold:
            continue
        if subtract_expect:
            arr = subtract_arr_expect(arr, c)
            arr = relu(arr)
        arr = arr[~np.isnan(arr)]
        entropy = scipy.stats.entropy(arr)
        bin_region = chrom, s, e
        grange = selector.binid_region2genome_range(bin_region)
        result = BedGraph(*grange, entropy)
        yield result


def write_bedgraph(bg_iter, outpath, mode='w'):
    """
    write the bedgraph to a file,
    bedgraph format, 4 columns:
        <chromosome>    <start>    <end>    <value>
    """
    with open(outpath, mode) as f:
        for chrom, start, end, value in bg_iter:
            fields = [str(i) for i in (chrom, start, end, value)]
            outline = "\t".join(fields) + "\n"
            f.write(outline)


def filter_abnormal(bg_iter):
    """
    filter out the abnormal bedgraphs
    """
    for bg in bg_iter:
        if bg.value == float('inf') or bg.value == -float('inf'):
            continue
        if not bg.chrom.startswith('chr'):
            bg = BedGraph("chr"+bg.chrom, bg.start, bg.end, bg.value)
        yield bg


def process_loci_chunk(selector, chrom, chunk, inner_window, outer_window, subtract_expect, coverage):
    results = sliding_cal_entropy(selector, chrom, chunk, inner_window, outer_window, subtract_expect, non_nan_threshold=coverage)
    filtered = filter_abnormal(results)
    return list(filtered)


def sort_bedGraph(path):
    """
    sort bedGraph files, yield sorted lines.
    """
    cmd = ['sort', '-k1,1', '-k2,2n', '-u', path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for line in p.stdout:
        yield line.decode('utf8')


def parse_bedgraph(line_iter):
    """
    yield sorted bedGraph tuples:
        (chrom, start, end, value)
    """
    for line in line_iter:
        chrom, start, end, value = line.strip().split()
        start, end = int(start), int(end)
        value = float(value)
        yield BedGraph(chrom, start, end, value)


def eliminate_overlap(bg_iter, window_size, overlap, resolution):
    """
    eliminate the overlap between blocks, use later value cover previous value.

        ----                        ----
      ----              ->        --
    ----                        --
    |--------------|            |--------------|
    """
    prev = next(bg_iter)
    for bg in bg_iter:
        if bg.chrom == prev.chrom:
            step = bg.start - prev.start
            if step <= 0:
                raise IOError("input bedgraph should be sorted.")
            else:
                if bg.start < prev.end:  # overlap
                    new_prev = BedGraph(prev.chrom, prev.start, bg.start, prev.value)
                    yield new_prev
                else:
                    yield prev
        else:
            yield prev
        prev = bg
    yield prev


@click.group()
def call_entropy():
    pass


@call_entropy.command()
@click.argument("cool-uri")
@click.argument("output")
@click.option("--window-size", "-w",
    type=int, default=1,
    show_default=True,
    help="The size of sliding window. "
         "The unit is number of bins.")
@click.option("--overlap", "-o",
    type=int, default=0,
    show_default=True,
    help="The size of sliding window overlap. "
         "Value should less than the window size, "
         "the unit is also the number of bins.")
@click.option("--inner-window", "-i",
    type=int, default=3,
    show_default=True,
    help="The size of inner window, "
         "interactions within inner window will be removed "
         "as self-ligation interaction. "
         "If not remove interactions within inner window, use 0. "
         "The unit is number of bins.")
@click.option("--outer-window", "-x",
    type=int, default=1001,
    show_default=True,
    help="The size of outer window, "
         "only consider the interactions within outer window. "
         "If consider interactions with whole chromosomes, use 0. "
         "The unit is number of bins.")
@click.option("--balance/--no-balance",
    default=True,
    show_default=True,
    help="Use balanced matrix or not.")
@click.option("--subtract-expect/--no-subtract-expect",
    default=True,
    show_default=True,
    help="Estimate the expect interaction value and minus it.")
@click.option("--coverage", "-c",
    default=0.5,
    show_default=True,
    help="The coverage rate threshold, "
         "only bins coverage large equal than this will be keeped. ")
@click.option("--processes", "-p",
    type=int, default=1,
    show_default=True,
    help="Number of processes.")
@click.option("--chunk-size", "-s",
    type=int, default=5000,
    show_default=True,
    help="How many blocks in one chunk, for parallel processing.")
def loci(cool_uri, output, window_size, overlap, inner_window, outer_window, balance, subtract_expect, coverage, processes, chunk_size):
    """

    \b
                         inner window
                         |-|
    |--------------------------------------------------------|
         |-                                 -|
                outer window
    \b
    Args
    ----
    cool_uri : str
        URI of input cool container.
    output : str
        Path to output BEDGRAPH file.
    """
    c = Cooler(cool_uri)
    resolution = c.info['bin-size']
    chromsizes = c.chromsizes.to_dict()
    chromsizes = {chr_ : (size//resolution) + 1 for chr_, size in chromsizes.items()}

    chr_chunks = chromosome_chunks(chromsizes, window_size, overlap, chunk_size)
    it1, it2 = tee(chr_chunks)  # split chr_chunks to chroms and chunks
    chroms = map(operator.itemgetter(0), it1)
    chunks = map(operator.itemgetter(1), it2)

    matrix_selector = MatrixSelector(c, balance=balance)

    if os.path.exists(output):
        subprocess.check_call(['rm', output])

    with ProcessPoolExecutor(max_workers=processes) as excuter:
        map_ = map if processes == 1 else excuter.map
        tmp_file = output + ".tmp"
        idx = 0
        args = repeat(matrix_selector), chroms, chunks, repeat(inner_window), repeat(outer_window), repeat(subtract_expect), repeat(coverage)
        for out_chunk in map_(process_loci_chunk, *args):
            print("chunk {}: {}/{} blocks".format(idx, len(out_chunk), chunk_size))
            write_bedgraph(out_chunk, tmp_file, mode='a')
            idx += 1

    sorted_lines = sort_bedGraph(tmp_file)
    bgs = parse_bedgraph(sorted_lines)
    non_overlap = eliminate_overlap(bgs, window_size, overlap, resolution)
    write_bedgraph(non_overlap, output)
    subprocess.check_call(['rm', output + ".tmp"])  # rm merged tmp file


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


def process_region_chunk(chunk, matrix_selector, non_nan_threshold=0.6):
    out_chunk = []
    for grange in chunk:
        m = matrix_selector.fetch(str(grange))
        arr = m[~np.isnan(m)]
        non_nan_rate = arr.shape[0] / (m.shape[0] * m.shape[1])
        if non_nan_rate < non_nan_threshold:
            continue
        val = scipy.stats.entropy(arr)
        out_chunk.append(BedGraph(*grange, val))
    return out_chunk


@call_entropy.command()
@click.argument("cool-uri")
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
@click.option("--coverage", "-c",
    default=0.5,
    show_default=True,
    help="The coverage rate threshold, "
         "only bins coverage large equal than this will be keeped. ")
@click.option("--processes", "-p",
    type=int, default=1,
    show_default=True,
    help="Number of processes.")
@click.option("--chunk-size", "-s",
    type=int, default=5000,
    show_default=True,
    help="How many blocks in one chunk, for parallel processing.")
def region(cool_uri, bed_path, output, inner_window, balance, coverage, processes, chunk_size):
    """
    \b
    Args
    ----
    cool_uri : str
        URI of input cool container.
    bed_path : str
        Path to input BED.
    output : str
        Path to output BEDGRAPH file.
    """
    c = Cooler(cool_uri)
    matrix_selector = MatrixSelector(c, balance=balance)
    regions = read_bed(bed_path)
    chunks = chunking(regions, chunk_size)
    if os.path.exists(output):
        subprocess.check_call(['rm', output])
    with ProcessPoolExecutor(max_workers=processes) as excuter:
        map_ = map if processes == 1 else excuter.map
        for out_chunk in map_(process_region_chunk, chunks, repeat(matrix_selector), repeat(coverage)):
            bgs = filter_abnormal(out_chunk)
            write_bedgraph(bgs, output, mode='a')


def unit_tests():
    def test_tiling():
        n = 10; w = 3; o = 2
        for b in tiling(n, w, o):
            print(b)

    def test_chunking():
        n = 10; w = 3; o = 2
        nc = 3
        blocks = tiling(n, w, o)
        chunks = chunking(blocks, nc)
        for c in chunks:
            print(c)

    def test_chromosome_chunks():
        chromsizes = {
            "1": 10,
            "2": 15,
        }
        w = 3; o = 2
        nc = 3
        chr_chunks = chromosome_chunks(chromsizes, w, o, nc)
        for chr_, chunk in chr_chunks:
            print(">" + chr_)
            for b in chunk:
                print(b)

    def test_eliminate_overlap():
        bgs = [
            BedGraph('1',    0, 1000, 1.0),
            BedGraph('1', 2000, 3000, 1.0),
            BedGraph('1', 2500, 3500, 1.1),
            BedGraph('1', 3000, 4000, 1.2),
            BedGraph('2', 1000, 2000, 1.0),
            BedGraph('2', 1500, 2500, 1.1),
            BedGraph('2', 2000, 3000, 1.2),
        ]
        res = eliminate_overlap((i for i in bgs), 2, 1, 500)
        for i in res:
            print(i)

    #test_tiling()
    #test_chunking()
    #test_chromosome_chunks()
    test_eliminate_overlap()


if __name__ == "__main__":
    eval("call_entropy()")
    #unit_tests()