# -*- coding: utf-8 -*-

"""
fastqreader
~~~~~~~~~~~
This module for read sequences from fastq files.

"""

from collections import Iterator
from itertools import islice


class FastqReader(Iterator):
    """
    Fastq file reader.
    read sequences from some fastq files.

    """
    def __init__(self, files, seqfilter=None):
        """
        files: file list of fastq file
        seqfilter: filter function,
            filter out the seq when seqfilter(seq) == False

        """
        if type(files) is list or type(files) is tuple:
            self.files = [open(f, 'r') for f in files]
        else:
            self.files = [open(files, 'r')]
        # point to the file which should be read
        self.pointer = 0
        self.current_file = self.files[0]
        # status repersent if there are files can be read or not.
        # 'U' means there are some can read, 'C' means all file readed
        self.status = 'U'
        self.filter = seqfilter

    def __del__(self):
        """ close the file handle """
        for f in self.files: 
            f.close()
            del f

    def next(self):
        """ behave as a iterator """
        lines = [self.current_file.readline() for i in range(4)]
        # arrived at the end of the file
        if lines == ['', '', '', '']:
            self.pointer += 1
            try:
                next_f = self.files[self.pointer]
            # all file readed
            except IndexError:
                self.status = 'C'
                raise StopIteration
            self.current_file = next_f
            return next(self)
        else:
            seq = lines[1].strip()
            # filter the seq by seqfilter
            if self.filter is None:
                return seq
            else:
                if self.filter(seq):
                    return seq
                else:
                    return next(self)

    def get_seqs(self, n):
        """ return a list of seqs. """
        return list(islice(self, n))
