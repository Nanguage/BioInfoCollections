# -*- coding: utf-8 -*-

"""
gtfparser
~~~~~~~~~
parse gtf file

useage:
>>> from gtfparser import GtfParser
>>> p = GtfParser("example.gtf")

# get item's information
>>> for item in p:
...     print(item.position)
...     print(item.gene_id)
...     print(item.strand)

"""

from collections import Iterator


class GtfParser(Iterator):
    def __init__(self, file):
        self.handle = open(file)

    def __del__(self):
        self.handle.close()

    def next(self):
        line = self.handle.readline()
        if line == '':
            raise StopIteration
        elif line.startswith('#'):
            return next(self)
        else:
            return GtfItem(line)
        

def _get_args(string):
    others = string.split(";")
    args = {}
    for i in others:
        i = i.strip()
        try:
            i = i.split(' ')
            args[i[0]] = i[1].strip('"')
        except IndexError:
            continue
    return args


class GtfItem:
    def __init__(self, line):
        self.line = line.strip()
        items = line.split("\t")
        self.source = items[1]
        self.type   = items[2]
        self.position = int(items[3]), int(items[4])
        self.strand = items[6]
        self.args = _get_args(items[-1])
        self.gene_id      = self.args['gene_id']
        
