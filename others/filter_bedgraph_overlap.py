#!/usr/bin/env python

"""
Filter out the overlaped bins in the bedGraph file.

usage:
    $ cat a.bedGraph | ./filter_bedgraph_overlap.py > a.filtered.bedGraph

"""

import fileinput


def print_line(items):
    print("\t".join([str(i) for i in items]))

def overlap(previous, current):
    """
    s_     e_
    |      |
        |     |
        s     e
    """
    chr, s, e, v = current
    chr_, s_, e_, v_ = previous
    if chr != chr_:
        return False
    if  s < e_:
        return True
    return False

previous = None

for line in fileinput.input():
    chr_, s, e, v = line.strip().split()
    s, e = int(s), int(e)
    items = chr_, s, e, v
    if not previous:
        print_line(items)
    else:
        if not overlap(previous, items):
            print_line(items)
    previous = items

