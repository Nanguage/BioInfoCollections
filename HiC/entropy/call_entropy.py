"""
"""

from concurrent.futures import ProcessPoolExecutor

import click
from cooler.api import Cooler
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



@click.command()
@click.argument("cool-uri")
@click.option("--window-size",
    type=int, default=1,
    help="The size of sliding window. " \
         "The unit is number of bins, " \
         "default 1")
@click.option("--overlap",
    type=int, default=0,
    help="The size of sliding window overlap. "\
         "Value should less than the window size, " \
         "the unit is also the number of bins. " \
         "default 0")
@click.option("--processes", "-n",
    type=int, default=1,
    help="Number of processes.")
@click.option("--chunk-size", "-c",
    type=int, default=10000)
def call_entropy(cool_uri, window_size, overlap, processes, chunk_size):
    c = Cooler(cool_uri)


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

    #test_tiling()
    test_chunking()


if __name__ == "__main__":
    #eval("call_entropy()")
    unit_tests()