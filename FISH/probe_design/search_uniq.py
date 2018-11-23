import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat, tee

import click
from Bio.SeqIO import FastaIO
from Bio.Blast.Applications import NcbiblastnCommandline


def save_fasta(seq_iter, path):
    with open(path, 'w') as f:
        for seq in seq_iter:
            seqname = seq.name
            if hasattr(seq, 'sub_range'):
                seqname += "|" + "{}-{}".format(*seq.sub_range)
            f.write(">" + seqname + "\n")
            f.write(str(seq.seq) + "\n")


def candidate_probes(seq, probe_length, search_step):
    """
    yield a series of candidate probe(sub-sequence)
    """
    start = 0
    end = start + probe_length
    if end > len(seq):
        seq.sub_range = (start, len(seq))
        yield seq
    else:
        while end <= len(seq):
            probe = seq[start:end]
            if len(probe) < probe_length:
                break
            probe.sub_range = (start, end)
            yield probe
            start += search_step
            end = start + probe_length


def is_unique_mapped(blastn_out_file, outfmt=6):
    """
    judge a blastn result is unique mapped or not.
    input a path to blastn result file in outfmt 6.
    """
    with open(blastn_out_file) as f:
        n = len(list(f))
    return n == 1


def passed(candidate_seq, blastn_db, tmpdir):
    """
    judge a candidate probe can pass(is unique mapped) or not.

    Parameters
    ----------
    candidate_seq : `Seq`
        candidate probe seq.
    blastn_db : str
        Path to blastn database.
    tmpdir : str
        Path to temporary dir, used for store the blasn result file.
    """
    seq = candidate_seq
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    seqname = seq.name+"_"+"{}-{}".format(*seq.sub_range)
    tmp_in  = os.path.join(tmpdir, seqname+".fa")
    save_fasta([seq], tmp_in)
    tmp_out = os.path.join(tmpdir, seqname+".tsv")
    cline = NcbiblastnCommandline(query=tmp_in, db=blastn_db, evalue=0.001, outfmt=6, out=tmp_out)
    stdout, stderr = cline()
    return is_unique_mapped(tmp_out)


def search_passed_probes(input_seqs, blastn_db, probe_length, search_step, blastn_tmpdir, threads):
    for seq in input_seqs:
        candidates = candidate_probes(seq, probe_length, search_step)
        candi_1, candi_2 = tee(candidates)
        with ProcessPoolExecutor(max_workers=threads) as executor:
            passed_iter = executor.map(passed, candi_2, repeat(blastn_db), repeat(blastn_tmpdir))
            for candidate, is_passed in zip(candi_1, passed_iter):
                print(candidate.name, candidate.sub_range)
                if is_passed:
                    yield candidate
                    break
            else:
                # not found passed
                print("[Warning] Not found passed probe in seq {}".format(seq.name))


@click.command()
@click.argument("input")
@click.argument("blastndb")
@click.argument("output")
@click.option("--probe-length", "-l",
    default=40,
    help="The length of probe(sub-sequence). default 40")
@click.option("--search-step", "-s",
    default=10,
    help="Unique map search step. default 10.")
@click.option("--blastn-tmpdir", "-b",
    default="./.blastn/",
    help="Temporary blastn directory. default .blastn")
@click.option("--threads", "-t",
    default=1,
    help="Use how many threads. default 1")
def main(input, blastndb, output, probe_length, search_step, blastn_tmpdir, threads):
    """
    Search unique mapped probe(sub-sequence)
    within a series of sequences stored in a fasta file.

    \b
    For example:
    select 30 candidate probe regions with length 500bp, firstly,
    ```
    $ python uniformly_spaced.py data/hg19.fa ./candidate.fa chr1:89000000-90000000 -n 30 -l 500
    ```
    then select unique maped probe(sub-sequence) from it.
    ```
    $ python search_uniq.py candidate.fa example/blastn_db/hg19 probe.fa
    ```

    \b
    Args
    ----
    input : str
        Path to input fasta file.
    blastndb : str
        Path to blastn database.
        build with `makeblastdb` command.
    output : str
        Path to output fasta file.
    
    """
    with open(input) as f:
        input_seqs = FastaIO.FastaIterator(f)
        probes = search_passed_probes(input_seqs, blastndb, probe_length, search_step, blastn_tmpdir, threads)
        save_fasta(probes, output)


if __name__ == "__main__":
    eval("main()")

