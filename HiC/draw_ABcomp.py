# -*- coding: utf-8 -*-

import sys
import argparse
from os.path import basename

import numpy as np
from scipy.interpolate import spline
import matplotlib.pyplot as plt


def argumrnt_parser():
    parser = argparse.ArgumentParser(
            description="draw A/B compartment hist(?) plots.")
    parser.add_argument("-o", help="output fig name")
    parser.add_argument("-i",
            help="input files, comma splited list, e.g. THP.txt,THPP.txt,Ra.txt")
    return parser


def parse_hist_file(fh):
    """
    parse hist file.
    The format of hist file:

    <Chromosome>    <start>    <end>    <count> ...

    """
    with fh as f:
        for line in f:
            chr_, s, e, c = line.strip().split()[:4]
            s, e = int(s), int(e)
            c = float(c)
            yield [chr_, s, e, c]


def list_avg_select(list_, n=10):
    length = len(list_) 
    span = length // n
    res = []
    for i in range(n):
        res.append(list_[i*span])
    res.append(list_[-1])
    return res


def plot(items, log=False, ax=plt.gca(), smooth=True, label=""):
    """
    plot bar graph
    """
    chrs, list_s, list_e, list_c = zip(*items)
    counts = np.asarray(list_c)
    num_bin = len(counts)
    ind = np.arange(num_bin) # x positions
    if log: 
        counts = np.log10(counts) # transform to log scale
    if smooth:
        newx = np.linspace(ind.min(), ind.max(), 100*len(ind))
        counts = spline(ind, counts, newx)
        ind = newx


    # fill between line and x axis
    x = ind
    y = counts
    ax.plot(ind, counts, linewidth=0)
    ax.fill_between(x, 0, y, where=y >= 0, facecolor='#ff9d9d', interpolate=True)
    ax.fill_between(x, 0, y, where=y <= 0, facecolor='#66ccff', interpolate=True)

    genome_len = list_e[-1]
    #plt.xticks([0, num_bin-1], [str(0), str(genome_len)])
    ind = list_avg_select(ind, 5)
    list_s = list_avg_select(list_s, 5)
    plt.xticks(ind, [str(i) for i in list_s])

    # remove right and top border
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # add label
    ax.text(0.5*x.min(), y.max(), label, fontsize=10, bbox={'facecolor':'white'})

    return ax
    

if __name__ == '__main__':
    parser = argumrnt_parser()
    args = parser.parse_args()

    files = args.i.split(",")
    labels = map(lambda s: basename(s).split('.')[0], files)

    is_log = False
    fig_name = args.o
    figsize = (20, 5)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1,
            figsize=figsize, sharex=True)

    fhs = map(open, files)
    for ax, fh, label in zip((ax1, ax2, ax3), fhs, labels):
        items = parse_hist_file(fh)
        plot(items, is_log, ax=ax, label=label)
    map(lambda fh:fh.close(), fhs)

    if fig_name != None:
        plt.savefig(fig_name)
    else:
        plt.show()
