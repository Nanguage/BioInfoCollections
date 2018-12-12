import re
from itertools import tee
from os.path import join
import os
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat, tee
from datetime import datetime
from collections import namedtuple

import numpy as np
import click
import pandas as pd
import h5py
from tqdm import tqdm
from cooler.api import Cooler


BED_FIELDS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb",
    "blockCount", "blockSizes", "blockStarts"
]


def read_bed(path):
    with open(path) as f:
        for line in f:
            bed_rec = line.strip().split()
            bed_rec[1] = int(bed_rec[1])
            bed_rec[2] = int(bed_rec[2])
            yield bed_rec


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
        resolution = self.cool.info['bin-size']
        start = binid1 * resolution
        end = binid2 * resolution
        grange = GenomeRange(chr_, start, end)
        return grange

    def genome_range2binid_region(self, genome_range):
        chr_, start, end = genome_range
        resolution = self.cool.info['bin-size']
        binid1 = start // resolution
        binid2 = end // resolution
        return (chr_, binid1, binid2)

    def confirm_chromosome_name(self, region):
        """
        confirm region's chromosome in cool

        Parameters
        ----------
        region : {str, `GenomeRange`}
        """
        if isinstance(region, str):
            grange = region_str2genome_range(region)
        else:
            grange = region
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


def count_range(mat_sel, chr, ref_pos, inner_window, up, down):
    """
    Count values in bigwig file within a genome range.

    Parameters
    ----------
    mat_sel : `MatrixSelector`
        matrix selector.
    chr : str
        chromosome name
    ref_pos : int 
        reference point position.
    inner_window : int
        Inner window size, unit: bin
    up : int
        how many bins up stream relative to reference point
    down : int
        down stream
    
    Return
    ------
    scores : `numpy.ndarray`
    """
    try:
        _, ref_pos, _ = mat_sel.genome_range2binid_region(GenomeRange(chr, ref_pos, ref_pos))
        outer_range = (chr, ref_pos - up, ref_pos + down + 1)
        flank = (inner_window-1) // 2
        inner_range = (chr, ref_pos - flank, ref_pos + flank + 1)
        arr = mat_sel.fetch_by_bin(inner_range, outer_range)
        arr = np.nanmean(arr, axis=0)
        scores = arr
    except (IOError, ValueError) as e:
#        print(e)
        arr = np.ones(shape=(up+down+1))
        arr[arr == 1] = np.nan
        scores = arr
    return scores


def iterator_to_dataframe(bed_iter):
    records = [i for i in bed_iter]
    columns = BED_FIELDS[:len(records[0])]
    df = pd.DataFrame(records, columns=columns)
    return df


def read_bed_df(path):
    df = pd.read_table(path, header=None, sep="\t")
    df.columns = BED_FIELDS[:len(df.columns)]
    return df


def split_uri(uri):
    if "::" not in uri:
        path = uri
        group = ""
    else:
        path, group = uri.split("::")
    return (path, group)


def dataframe_to_hdf5(df, base_uri, gname):
    path, group = split_uri(base_uri)
    hdf5_block_open(path).close()  # waiting for hdf lock release
    with pd.HDFStore(path) as store:
        group = join(group, gname)
        store[group] = df


def hdf5_block_open(path, wait_time=5):
    """
    open the hdf file.
    wait, if it is be locked.
    """
    while True:
        try:
            f = h5py.File(path)
            break
        except OSError as oe:
            print(str(oe))
            print("[Warning] {} is locked waiting for lock release.".format(path))
            time.sleep(wait_time)
    return f


def incremental_write_h5dset(base_uri, dset_name, array):
    path, group = split_uri(base_uri)
    wait_time = 5
    f = hdf5_block_open(path, wait_time)
    grp = f[group] if group else f
    if dset_name not in grp:
        grp.create_dataset(dset_name, data=array, maxshape=(None, array.shape[1]))
    else:
        dset = grp[dset_name]
        chunk_size = array.shape[0]
        dset.resize(dset.shape[0] + chunk_size, axis=0)
        dset[-chunk_size:] = array
    f.close()


def clean_dataset(base_uri, dset_name):
    path, group = split_uri(base_uri)
    with hdf5_block_open(path) as f:
        grp = f[group] if group else f
        if dset_name in grp:  # clean dataset
            del grp[dset_name]


def scores_iter_to_hdf5(scores_iter, base_uri, dset_name, chunk_size=2000, dtype='float64'):
    clean_dataset(base_uri, dset_name)
    chunk = []
    for idx, scores in tqdm(enumerate(scores_iter), total=num_records):
        if (idx != 0) and (idx % chunk_size == 0):
            chunk_arr = np.asarray(chunk, dtype=dtype)
            incremental_write_h5dset(base_uri, dset_name, chunk_arr)
            chunk = []
        chunk.append(scores)
    if chunk:
        chunk_arr = np.asarray(chunk, dtype=dtype)
        incremental_write_h5dset(base_uri, dset_name, chunk_arr)


def write_meta_info(h5group_uri, bed, cool_uri, inner_window, up_stream, down_stream):
    path, group = split_uri(h5group_uri)
    f = hdf5_block_open(path)
    now = str(datetime.now())
    grp = f[group] if group else f
    grp.attrs.update({
        'create_date': now,
        'reference_bed': bed,
        'source_cool': cool_uri,
        'up_stream_bins': up_stream,
        'down_stream_bins': down_stream,
        'inner_window_size': inner_window,
    })
    f.close()


@click.command()
@click.argument("bed")
@click.argument("cool_uri")
@click.argument("h5group_uri")
@click.option("--inner-window", "-i",
    default=3,
    show_default=True,
    help="The inner window size, unit: bin.")
@click.option("--up-stream", "-u",
    default=1000,
    show_default=True,
    help="Up stream range, unit: bin")
@click.option("--down-stream", "-d",
    default=1000,
    show_default=True,
    help="Down stream range, unit: bin")
@click.option("--balance/--no-balance",
    default=True,
    show_default=True,
    help="Use balanced matrix or not.")
@click.option("--processes", "-p",
    default=1,
    show_default=True,
    help="How many process to run.")
def stats_v4c(bed, cool_uri, h5group_uri, inner_window, up_stream, down_stream, balance, processes):
    """
    Compute the value matrix in bigwig around start position in bed file.

    \b
    Args
    ----
    bed : str
        Path to input bed file.
    cool_uri : str
        URI to cool.
    h5group_uri : str
        URI of output HDF5 file group, like:
            ./test.h5
            ./test.h5::/virtual4c/
    """
    path, group = split_uri(h5group_uri)
    if not os.path.exists(path):
        h5py.File(path).close()  # create, if file not exist.
    df = read_bed_df(bed)
    global num_records
    num_records = df.shape[0]  # for create progress bar
    dataframe_to_hdf5(df, h5group_uri, "ref_bed")
    bed_recs = read_bed(bed)
    cool = Cooler(cool_uri)
    mat_sel = MatrixSelector(cool, balance=balance)
    def iterover_fetch_scores(iter):
        chrs, ref_pos = tee(iter)
        chrs = (rec[0] for rec in chrs)
        ref_pos = (rec[1] for rec in ref_pos)
        map_ = ProcessPoolExecutor(max_workers=processes).map if processes > 1 else map
        args = (repeat(mat_sel), chrs, ref_pos, repeat(inner_window), repeat(up_stream), repeat(down_stream))
        for scores in map_(count_range, *args):
            yield scores
    scores_iter = iterover_fetch_scores(bed_recs)
    incremental_chunk_size = 20
    scores_iter_to_hdf5(scores_iter, h5group_uri, "matrix", incremental_chunk_size)
    write_meta_info(h5group_uri, bed, cool_uri, inner_window, up_stream, down_stream)


if __name__ == "__main__":
    eval("stats_v4c()")
