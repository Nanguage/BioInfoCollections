from itertools import tee
from os.path import join
import os

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


def count_range(bw, chr, ref_pos, up, down, n_bins):
    chroms = bw.chroms()
    if chr not in chroms:
        chr = change_chr_name(chr)
        if chr not in chroms:
            return
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
    with pd.HDFStore(path) as store:
        group = join(group, gname)
        store[group] = df


def incremental_write_h5dset(group_handler, dset_name, array):
    grp = group_handler
    if dset_name not in grp:
        grp.create_dataset(dset_name, data=array, maxshape=(None, array.shape[1]))
    else:
        dset = grp[dset_name]
        chunk_size = array.shape[0]
        dset.resize(dset.shape[0] + chunk_size, axis=0)
        dset[-chunk_size:] = array


def arr_iter_to_hdf5(arr_iter, base_uri, dset_name, chunk_size=2000, dtype='float64'):
    path, group = split_uri(base_uri)
    with h5py.File(path) as f:
        grp = f[group] if group else f
        if dset_name in grp:  # clean dataset
            del grp[dset_name]
        chunk = []
        for idx, arr in tqdm(enumerate(arr_iter), total=num_records):
            if (idx != 0) and (idx % chunk_size == 0):
                chunk_arr = np.asarray(chunk, dtype=dtype)
#                import ipdb; ipdb.set_trace()
                incremental_write_h5dset(grp, dset_name, chunk_arr)
                chunk = []
            chunk.append(arr)
        if chunk:
            chunk_arr = np.asarray(chunk, dtype=dtype)
            incremental_write_h5dset(grp, dset_name, chunk_arr)


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
def compute_matrix(bed, bigwig, h5group_uri, num_bins, up_stream, down_stream):
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
    bw = pyBigWig.open(bigwig)
    bed_recs = read_bed(bed)
    def iterover_fetch_array(iter):
        for rec in iter:
            arr = count_range(bw, rec[0], rec[1], up_stream, down_stream, num_bins)
            yield arr
    arr_iter = iterover_fetch_array(bed_recs)
    incremental_chunk_size = 20
    arr_iter_to_hdf5(arr_iter, h5group_uri, "matrix", incremental_chunk_size)


if __name__ == "__main__":
    eval("stats_bigwig()")
