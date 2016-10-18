# -*- coding: utf-8 -*-

"""
fastqreader
~~~~~~
This module for read sequences from fastq files.

"""

from collections import Iterator

class FastqReader(Iterator):
    '''
    Fastq file reader.
    read sequences from some fastq files.

    '''
    def __init__(self, files, seqfilter=None):
        '''
        files: file list of fastq file
        seqfilter: filter function,
            filter out the seq when seqfilter(seq) == False

        '''
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
        '''close the file handle'''
        for f in self.files: 
            f.close()
            del f

    def get_seqs(self, step=1000000):
        '''
        read some sequence 
        step: how many lines once readed
        
        '''
        result_seqs = []
        for i in xrange(step):
            # arrived at the end of the file, and jump over one line
            if self.current_file.readline() == '':
                self.pointer += 1
                try:
                    next_f = self.files[self.pointer]
                # all file readed
                except IndexError:
                    self.status = 'C'
                    return result_seqs
                self.current_file = next_f
            seq = self.current_file.readline().strip()
            # filter the seq by seqfilter
            if self.filter is not None:
                if self.filter(seq):
                    result_seqs.append(seq)
            else:
                result_seqs.append(seq)
            # jump over 2 lines which is not seq
            for i in xrange(2): self.current_file.readline()
        return result_seqs

    def get_seq(self):
        '''read one sequence'''
        return self.get_seqs(step=1)

    def next(self):
        ''' behave as a iterator '''
        result = None
        lines = [self.current_file.readline() for i in range(4)]
        # arrived at the end of the file
        if lines == '':
            self.pointer += 1
            try:
                next_f = self.files[self.pointer]
            # all file readed
            except IndexError:
                self.status = 'C'
                raise StopIteration
            self.current_file = next_f
        seq = lines[1]
        # filter the seq by seqfilter
        if self.filter is not None:
            if self.filter(seq):
                result = seq
        else:
            result = seq
        return result


