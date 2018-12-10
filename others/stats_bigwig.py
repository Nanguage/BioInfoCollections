from itertools import tee
from os.path import join
import os
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat, tee

import numpy as np
import click
import pyBigWig
import pandas as pd
import h5py
from tqdm import tqdm


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


def write_bed(bed_iter, path):
    with open(path, 'w') as f:
        for rec in bed_iter:
            outline = "\t".join([str(i) for i in rec])
            f.write(outline + "\n")


def change_chr_name(chr_name):
    if chr_name.startswith("chr"):
        return chr_name.replace("chr", "")
    else:
        return "chr" + chr_name


def count_range(bigwig, chr, ref_pos, up, down, n_bins):
    """
    Count values in bigwig file within a genome range.

    Parameters
    ----------
    bigwig : str
        Path to bigwig file.
    chr : str
        chromosome name
    ref_pos : int 
        reference point position.
    up : int
        how many bp up stream relative to reference point
    down : int
        down stream
    n_bins: int
        how many bins.
    
    Return
    ------
    scores : list of {None, float}
    """
    bw = pyBigWig.open(bigwig)
    chroms = bw.chroms()
    if chr not in chroms:
        chr = change_chr_name(chr)
        if chr not in chroms:
            return [None] * n_bins
    s_, e_ = ref_pos - up, ref_pos + down
    chr_len = chroms[chr]
    x1 = max(0, s_)
    x2 = min(chr_len, e_)
    scores = bw.stats(chr, x1, x2, nBins=n_bins)
    if len(scores) != n_bins:
        len_per_bin = (up + down) / n_bins
        if s_ < 0:
            n_pre = int(abs(s_) / len_per_bin)
            scores = [None] * n_pre + scores
        if e_ > chr_len:
            n_post = int((e_ - chr_len) / len_per_bin)
            scores = scores + [None] * n_post 
        assert len(scores) == n_bins
    return scores


def iterator_to_dataframe(bed_iter):
    records = [i for i in bed_iter]
    columns = BED_FIELDS[:len(records[0])]
    df = pd.DataFrame(records, columns=columns)
    return df


@click.group()
def stats_bigwig():
    pass


@stats_bigwig.command()
@click.argument("bed")
@click.argument("bigwig")
@click.argument("output")
@click.option("--up-stream", "-u",
    default=1000,
    show_default=True,
    help="Up stream range")
@click.option("--down-stream", "-d",
    default=1000,
    show_default=True,
    help="Down stream range")
def region_mean(bed, bigwig, output, up_stream, down_stream):
    """

    Count the mean value in bigwig file, around start position in bed file.

    \b
    Args
    ----
    bed : str
        Path to input bed file.
    bigwig : str
        Path to input bigwig file.
    output : str
        Path to output bed file.
    """
    bw = pyBigWig.open(bigwig)
    def process_bed(iter):
        for rec in iter:
            ref_pos = rec[1]
            score = count_range(bw, rec[0], ref_pos, up_stream, down_stream, n_bins=1)[0]
            score = score or '.'
            rec[4] = score
            yield rec
    bed_recs = read_bed(bed)
    scored_bed = process_bed(bed_recs)
    write_bed(scored_bed, output)


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
    with h5py.File(path) as f:
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


@stats_bigwig.command()
@click.argument("bed")
@click.argument("bigwig")
@click.argument("h5group_uri")
@click.option("--num_bins", "-n",
    default=1000,
    help="Number bins within one fetched region.")
@click.option("--up-stream", "-u",
    default=1000,
    show_default=True,
    help="Up stream range")
@click.option("--down-stream", "-d",
    default=1000,
    show_default=True,
    help="Down stream range")
@click.option("--processes", "-p",
    default=1,
    show_default=True,
    help="How many process to run.")
def compute_matrix(bed, bigwig, h5group_uri, num_bins, up_stream, down_stream, processes):
    """
    Compute the value matrix in bigwig around start position in bed file.

    \b
    Args
    ----
    bed : str
        Path to input bed file.
    bigwig : str
        Path to input bigwig file.
    h5group_uri : str
        URI of output HDF5 file group, like:
            ./test.h5
            ./test.h5::/h3k27ac/
    """
    path, group = split_uri(h5group_uri)
    if not os.path.exists(path):
        f = h5py.File(path)
        f.close()
    df = read_bed_df(bed)
    global num_records
    num_records = df.shape[0]
    dataframe_to_hdf5(df, h5group_uri, "ref_bed")
    bed_recs = read_bed(bed)
    def iterover_fetch_scores(iter):
        chrs, ref_pos = tee(iter)
        chrs = (rec[0] for rec in chrs)
        ref_pos = (rec[1] for rec in ref_pos)
        map_ = ProcessPoolExecutor(max_workers=processes).map if processes > 1 else map
        args = (repeat(bigwig), chrs, ref_pos, repeat(up_stream), repeat(down_stream), repeat(num_bins))
        for scores in map_(count_range, *args):
            yield scores
    scores_iter = iterover_fetch_scores(bed_recs)
    incremental_chunk_size = 20
    scores_iter_to_hdf5(scores_iter, h5group_uri, "matrix", incremental_chunk_size)


if __name__ == "__main__":
    eval("stats_bigwig()")
