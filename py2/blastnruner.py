# -*- coding: utf-8 -*-

"""
blastnruner
~~~~~~~~~~~
This module for do Blast(blastn) and return the match information.

"""

import os
import random

class BlastnRuner():
    """
    Do blastn and parse the result, return the matches.

    """
    def __init__(self, blastdb, outfmt=6, evalue=0.01, short=False, cache='./blastncache'):
        """
        blastdb: path to the database of blast
        outfmt: output format of blastn, default 5 is xml
        evalue: describes the number of hits one can "expect" to see
            by chance when searching a database of a particular size.
        short: if True, will add '-task blastn-short to command'
            use when sequence length less than 30 nucleotides.
        cache: cache directory.

        """
        self.id = str(random.randint(10000, 20000)) + str(id(self)) + str(os.getpid())
        self._fmt = outfmt
        # buffer for storage sequences
        self.buffer  = []
        # make it if cache directory not exist
        if not os.path.exists(cache):
            os.mkdir(cache)
        self.cache = cache
        # in and out file's name of blastn
        self.in_file  = os.path.join(cache, self.id + '.in')
        self.out_file = os.path.join(cache, self.id + '.out')
        # blastn command templament
        self.command = "blastn -query {{in_file}} -db {db} -out {{out_file}} -outfmt {fmt} -evalue {e}".format(db=blastdb, fmt=outfmt, e=evalue)
        if short: self.command += ' -task blastn-short'

    def __del__(self):
        """Clean the cache file."""
        if os.path.exists(self.in_file):
            os.remove(self.in_file)
        if os.path.exists(self.out_file):
            os.remove(self.out_file)

    def add(self, seqs):
        """Add sequences to buffer."""
        if type(seqs) is list:
            self.buffer += seqs
        else:
            self.buffer.append(seqs)

    def _write(self):
        """Write the sequences in the buffer to in_file."""
        with open(self.in_file, 'w') as f:
            for i, s in enumerate(self.buffer):
                f.write('>%d\n'%i)
                f.write('%s\n'%s)

    def _read(self):
        """Read blastn output."""
        with open(self.out_file, 'r') as f:
            result = f.read()
        return result

    def parse_xml(self, xml):
        """
        Parse the output xml.
        xml: struct(useful) like this:
            <Iteration>
                <Iteration_iter-num>1</Iteration_iter-num>
                <Iteration_query-len>12</Iteration_query-len>
                <Iteration_hits>
                    <Hit></Hit>
                </Iteration_hits>
            </Iteration>

        """
        # TODO
        raise NotImplementedError("Can only parse fmt 6 now")

    def parse_fmt6(self, fmt6):
        """
        Parse NO.6 output format.
        fmt6:
          format:
            Query Target Percent Alignment mismatch gaps start_q end_q start_t end_t e-value bit_score
          example:
            1 Rv1023  100.00  12  0   0   1   12  16  5   0.002   24.3
            3 Rv0911  100.00  12  0   0   1   12  14  3   0.002   24.3
            3 Rv2518c 100.00  11  0   0   2   12  13  3   0.010   22.3
            ...
        RETURN: a list of matches like:
            [
                [('Rv1023, 100, 12, ...)],
                None,
                [('Rv0911', 100, 12, ...), ('Rv2518c', 100, 11, ...)],
                ...
            ]
            one item corresponding to one sequence in buffer
            if one item can't find match the item is: None

        """
        # No matched result
        if fmt6 == '':
            return [None for i in self.buffer]
        result = []
        lines = fmt6.split('\n')
        lines = [line.split('\t') for line in lines]
        # make a dict mapping query id to it's other match infomation.
        query2info = {}
        for l in lines:
            query2info.setdefault(int(l[0]), [])
            query2info[int(l[0])].append(l[1:])
        # set result for every sequence in the buffer.
        for i, s in enumerate(self.buffer):
            info = query2info.get(i, None)
            result.append(info if info else None)
        return result

    def blastn(self):
        """
        The interface of the class.
        CALL: _write->_read->parse_
        RETURN: a list of best match

        """
        self._write()
        command = self.command.format(in_file=self.in_file,\
                out_file=self.out_file)
        os.system(command)
        bn_result = self._read().strip()
        if self._fmt == 6:
            result = self.parse_fmt6(bn_result)
        elif self._fmt == 5:
            result =  self.parse_xml(bn_result)
        # clean buffer
        self.buffer = []
        return result

    def blastn_one(self, seq):
        """
        Blast one sequence use subprocess.

        """
        # TODO
        raise NotImplementedError("unimplemented method.")
