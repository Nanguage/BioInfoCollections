import sys

import pandas as pd
import click


def read_tad_file(path):
    df = pd.read_table(path, sep="\t", header=None)
    columns = list(df.columns)
    columns[:3] = ['chr', 'start', 'end']
    df.columns = columns
    return df


@click.command()
@click.argument("tad1")
@click.argument("tad2")
@click.argument("output")
@click.option("--diff-threshold", "-d",
    default=10000,
    show_default=True,
    help="Threshold position diffenence.")
@click.option("--expand/--no-expand",
    default=True,
    help="Output similar tads in tad2 file.")
def call_conserved(tad1, tad2, output, diff_threshold, expand):
    """
    Find the conserved TADs between two sample.

    \b
    Parameters
    ----------
    TAD1 : str
        Path to TAD file of sample1
    TAD2 : str
        Path to TAD file of sample2
    output : str
        Path to output file, use '-' output to stdout.
    """
    span = diff_threshold
    df1 = read_tad_file(tad1)
    df2 = read_tad_file(tad2)
    conserved = []

    out_f = sys.stdout if output == "-" else open(output, 'w')
    
    for _, row in df1.iterrows():
        chr_, start, end = row.chr, row.start, row.end
        start_range = (start - span, start + span)
        end_range = (end - span, end + span)
        candidates = df2[(df2.chr == chr_) &
                         (df2.start >= start_range[0]) & (df2.start <= start_range[1]) &
                         (df2.end   >= end_range[0])   & (df2.end   <= end_range[1])]
        if candidates.shape[0] > 0:
            conserved.append(row)
            if expand:
                print(*list(row), sep="\t", end="\t", file=out_f)
                candi_ = [list(c) for _, c in candidates.iterrows()]
                print("\t".join([str(i) for subl in candi_ for i in subl]), sep="\t", file=out_f)
            else:
                print(*list(row), sep="\t", file=out_f)

    if not out_f is sys.stdout:
        out_f.close()

if __name__ == "__main__":
    eval("call_conserved()")
