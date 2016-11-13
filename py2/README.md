bioinfo collections (python2)
=============================

### [fastqreader](https://github.com/Nanguage/BioInfoCollections/tree/master/py2/fastqreader.py)
A simple fastq file format reader, provide a iterable and continuous reader for multiple fastq files. And you can add a sequence filter for it.
```
>>> from fastqreader import FastqReader
>>> r = FastqReader(['./sample1.fastq', './sample2.fastq', './sample3.fastq'])
>>> for seq in r: print seq
ATTTCCGGGTTTTTTCCCC
ACTTCCGAGTTTTTTCCCC
GCTTCCGAGTTTTTTCCCT
...
AATTCGGCCTTACCA
>>> r = FastqReader(['./sample1.fastq', './sample2.fastq', './sample3.fastq'], seqfilter=lambda s:len(s)<=10)
>>> for seq in r: print seq
ATTCGG
ATTAGC
AATC
...
AAAGGTC
```

### [blastnruner](https://github.com/Nanguage/BioInfoCollections/tree/master/py2/blastnruner.py)
The wrap of ncbi blastn, provide a convenient blastn interface to python.
```
>>> seqs = ['ATTTCCGGGTTTTTTCCCC', 'ACTTCCGAGTTTTTTCCCC', 'GCTTCCGAGTTTTTTCCCT']
>>> from blastnruner import BlastnRuner
>>> blaster = BlastnRuner()
>>> for s in seqs: 
...     blaster.add(s)
>>> res = blaster.blastn()
>>> for r in res:
...     print r

['Rv1024', '100.00', '12', '0', '0', '1', '12', '16', '5', '0.002', '24.3']
None
['Rv2518c', '100.00', '11', '0', '0', '2', '12', '13', '3', '0.010', '22.3']
```
