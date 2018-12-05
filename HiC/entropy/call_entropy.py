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

    def binid_region2genome_range(self, binid_region):
        chr_, binid1, binid2 = binid_region
        assert chr_ in self.cool.chromnames
        resolution = self.cool.info['bin-size']
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
    eliminate inner window of fetched array in loci entropy calling.
    """
    center = (start + end) // 2
    flank = (inner_window -1) // 2
    len_ = mat.shape[1]
    s = max(0, center - flank)
    e = min(center + flank+1, len_)
    mat[:, s:e] = np.nan
    return mat


def sliding_cal_entropy(selector, chrom, chunk, inner_window, non_nan_threshold=0.6):
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
    non_nan_threshold : float
        non-NaN value threshold, if less than this will be droped.
    """
    start, end = chunk[0][0], chunk[-1][-1]
    matrix = selector.fetch_by_bin((chrom, start, end), chrom)
    for s, e in chunk:
        s_, e_ = s - start, e - start
        m = matrix[s_:e_, :]
        if inner_window > 0:
            m = eliminate_inner(m, s, e, inner_window)
        arr = m[~np.isnan(m)]
        non_nan_rate = arr.shape[0] / (m.shape[0] * m.shape[1])
        if non_nan_rate < non_nan_threshold:
            continue
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


def process_loci_chunk(selector, chrom, chunk, output, inner_window, coverage):
    results = sliding_cal_entropy(selector, chrom, chunk, inner_window, non_nan_threshold=coverage)
    filtered = filter_abnormal(results)
    try_id = 0
    tmp_path = lambda id_: output + ".tmp.{}".format(id_)
    while os.path.exists(tmp_path(try_id)):
        try_id += 1
    tmp_file = tmp_path(try_id)
    write_bedgraph(filtered, tmp_file)
    return tmp_file


def merge_and_sort(tmp_files):
    """
    Merge tmp bedGraph files, and sort it, 
    yield sorted lines.
    """
    merged_file = re.sub("\.tmp\..", "", tmp_files[0]) + ".tmp"
    # merge tmp_files
    cmd = "cat " + " ".join(tmp_files) + " > " + merged_file
    subprocess.check_call(cmd, shell=True)
    # rm tmp files
    cmd = ['rm'] + tmp_files
    subprocess.check_call(cmd)
    # sort output
    cmd = ['sort', '-k1,1', '-k2,2n', '-u', merged_file]
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
            expected_step = (window_size - overlap) * resolution
            step = bg.start - prev.start
            if step <= 0:
                raise IOError("input bedgraph should be sorted.")
            else:
                if step != expected_step:
                    yield prev
                else:
                    new_prev = BedGraph(prev.chrom, prev.start, bg.start, prev.value)
                    yield new_prev
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
    type=int, default=10000,
    show_default=True,
    help="How many blocks in one chunk, for parallel processing.")
def loci(cool_uri, output, window_size, overlap, inner_window, balance, coverage, processes, chunk_size):
    """
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

    with ProcessPoolExecutor(max_workers=processes) as excuter:
        tmp_files = []
        map_ = map if processes == 1 else excuter.map
        for fname in map_(process_loci_chunk, repeat(matrix_selector), chroms, chunks, repeat(output), repeat(inner_window), repeat(coverage)):
            tmp_files.append(fname)

    sorted_lines = merge_and_sort(tmp_files)
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