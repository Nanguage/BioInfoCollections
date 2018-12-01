import pandas as pd
import click

def skip_lines(path):
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith("track"):
                n += 1
            else:
                break
    return n


def read_bed(path):
    n_skip = skip_lines(path)
    df = pd.read_table(path, sep="\t", header=None, skiprows=n_skip)
    base_cols = ['chr', 'start', 'end']
    n_col = len(df.columns)
    if n_col == 4:
        columns = base_cols + ['value']
    else:
        columns = base_cols
        if n_col >= 6:
            columns += ['name', 'score', 'strand']
        if n_col >= 9:
            columns += ['thickStart', 'thickEnd', 'itemRgb']
        if n_col == 12:
            columns += ['blockCount', 'blockSizes', 'blockStarts']
    df.columns = columns
    return df


def region_str(df):
    ser = df.chr + '_' + df.start.map(str) + '_' + df.end.map(str)
    return ser


@click.command()
@click.argument("bedgraph")
@click.argument("output")
@click.option("--ref-bed", "-r",
    help="reference BED file.")
def bedgraph2bed(bedgraph, output, ref_bed):
    """
    Expand bedGraph to BED.

    Default expand to BED6, if set reference BED, 
    substitude the value section with bedgraph value.
    """
    bg = read_bed(bedgraph)
    if ref_bed:
        ref_bed = read_bed(ref_bed)
        outbed = ref_bed
        bg.index = region_str(bg)
        outbed.index = region_str(outbed)
        outbed = outbed.loc[bg.index]
        outbed.score = bg.value
    else:
        outbed = bg
        outbed['name'] = '.'
        outbed['score'] = bg.value
        outbed['strand'] = '.'
        outbed = outbed[['chr', 'start', 'end', 'name', 'score', 'strand']]
    outbed.to_csv(output, header=False, sep="\t", index=False)


if __name__ == "__main__":
    eval("bedgraph2bed()")

