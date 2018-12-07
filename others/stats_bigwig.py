import click
import pyBigWig


def read_bed(path):
    with open(path) as f:
        for line in f:
            bed_rec = line.strip().split()
            bed_rec[1] = int(bed_rec[1])
            bed_rec[2] = int(bed_rec[2])
            yield bed_rec


def write_bed(rec_iter, path):
    with open(path, 'w') as f:
        for rec in rec_iter:
            outline = "\t".join([str(i) for i in rec])
            f.write(outline + "\n")


def change_chr_name(chr_name):
    if chr_name.startswith("chr"):
        return chr_name.replace("chr", "")
    else:
        return "chr" + chr_name


def count_range_mean(bw, chr, ref_pos, up, down):
    chroms = bw.chroms()
    if chr not in chroms:
        chr = change_chr_name(chr)
        if chr not in chroms:
            return
    x1 = max(0, ref_pos - up)
    x2 = min(chroms[chr], ref_pos + down)
    return bw.stats(chr, x1, x2)[0]


@click.command()
@click.argument("bed")
@click.argument("bigwig")
@click.argument("output")
@click.option("--up-stream", "-u",
    default=1000,
    show_default=True,
    help="Up stream range")
@click.option("--down-stream", "-d",
    default=1000,
    show_default=True,
    help="Down stream range")
def stats_bigwig(bed, bigwig, output, up_stream, down_stream):
    """

    Count the mean value in bigwig file, around start position in bed file.

    \b
    Args
    ----
    bed : str
        Path to input bed file.
    bigwig : str
        Path to input bigwig file.
    output : str
        Path to output bed file.
    """
    bw = pyBigWig.open(bigwig)
    def process_bed(iter):
        for rec in iter:
            ref_pos = rec[1]
            score = count_range_mean(bw, rec[0], ref_pos, up_stream, down_stream)
            score = score or '.'
            rec[4] = score
            yield rec
    bed_recs = read_bed(bed)
    scored_bed = process_bed(bed_recs)
    write_bed(scored_bed, output)


if __name__ == "__main__":
    eval("stats_bigwig()")
