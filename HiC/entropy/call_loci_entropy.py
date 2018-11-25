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


class MatrixSelector(object):
    """
    Selector for fetch the matrix from cool file.
    """
    def __init__(self, cool, balance=True):
        self.cool = cool
        self.balance = balance

    def binid_region2genome_region(self, binid_region):
        chr_, binid1, binid2 = binid_region
        assert chr_ in self.cool.chromnames
        resolution = self.cool.info['bin-size']
        start = binid1 * resolution
        end = binid2 * resolution
        region = (chr_, start, end)
        return region

    def fetch(self, bin_region_1, bin_region_2):
        m = self.cool.matrix(balance=self.balance)
        def format_region(chr_, start, end):
            return "{}:{}-{}".format(chr_, start, end)
        def convert_region(bin_region):
            if isinstance(bin_region, tuple):
                region = self.binid_region2genome_region(bin_region)
                region = format_region(*region)
            else:
                assert bin_region in self.cool.chromnames
                region = bin_region
            return region
        region1 = convert_region(bin_region_1)
        region2 = convert_region(bin_region_2)
        return m.fetch(region1, region2)


BedGraph = namedtuple("BedGraph", ['chrom', 'start', 'end', 'value'])


def sliding_cal_entropy(selector, chrom, chunk, non_nan_threshold=0.6):
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
    """
    start, end = chunk[0][0], chunk[-1][-1]
    matrix = selector.fetch((chrom, start, end), chrom)
    for s, e in chunk:
        s_, e_ = s - start, e - start
        m = matrix[s_:e_, :]
        arr = m[~np.isnan(m)]
        non_nan_rate = arr.shape[0] / (m.shape[0] * m.shape[1])
        if non_nan_rate < non_nan_threshold:
            continue
        entropy = scipy.stats.entropy(arr)
        bin_region = chrom, s, e
        region = selector.binid_region2genome_region(bin_region)
        result = BedGraph(*region, entropy)
        yield result


def write_bedgraph(bg_iter, outpath):
    """
    write the bedgraph to a file,
    bedgraph format, 4 columns:
        <chromosome>    <start>    <end>    <value>
    """
    with open(outpath, 'w') as f:
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


def process_chunk(selector, chrom, chunk, output, coverage):
    results = sliding_cal_entropy(selector, chrom, chunk, non_nan_threshold=coverage)
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


@click.command()
@click.argument("cool-uri")
@click.argument("output")
@click.option("--window-size", "-w",
    type=int, default=1,
    help="The size of sliding window. " \
         "The unit is number of bins, " \
         "default 1")
@click.option("--overlap", "-o",
    type=int, default=0,
    help="The size of sliding window overlap. "\
         "Value should less than the window size, " \
         "the unit is also the number of bins. " \
         "default 0")
@click.option("--balance/--no-balance",
    default=True,
    help="Use balanced matrix or not.")
@click.option("--processes", "-p",
    type=int, default=1,
    help="Number of processes.")
@click.option("--chunk-size", "-s",
    type=int, default=10000,
    help="How many blocks in one chunk, for parallel processing.")
@click.option("--coverage", "-c",
    default=0.5,
    help="The coverage rate threshold, " + \
         "only bins coverage large equal than this will be keeped. " + \
         "default 0.5")
def call_loci_entropy(cool_uri, output, window_size, overlap, balance, processes, chunk_size, coverage):
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
        for fname in excuter.map(process_chunk, repeat(matrix_selector), chroms, chunks, repeat(output), repeat(coverage)):
            tmp_files.append(fname)

    sorted_lines = merge_and_sort(tmp_files)
    bgs = parse_bedgraph(sorted_lines)
    non_overlap = eliminate_overlap(bgs, window_size, overlap, resolution)
    write_bedgraph(non_overlap, output)
    subprocess.check_call(['rm', output + ".tmp"])  # rm merged tmp file


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
    eval("call_loci_entropy()")
    #unit_tests()