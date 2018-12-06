import re
from collections import namedtuple

import click
import pandas as pd

GTF_FIELDS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes", "comments"]
BED6_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]


def infer_skip_header_rows(path):
    i = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                i += 1
            else:
                break
    return i


def read_gtf(path):
    skip = infer_skip_header_rows(path)
    gtf = pd.read_table(path, sep="\t", header=None, skiprows=skip)
    gtf.columns = GTF_FIELDS[:len(gtf.columns)]
    return gtf


def extract_section_from_attributes(section_name, attributes):
    return attributes.str.extract('.*{} "(.*?)"'.format(section_name), expand=False)


def gtfdf2beddf(gtf, name):
    name = extract_section_from_attributes(name, gtf.attributes)
    bed = pd.DataFrame({
        "chrom": gtf.seqname,
        "start": gtf.start,
        "end": gtf.end,
        "name": name,
        "score": gtf.score,
        "strand": gtf.strand,
    })
    bed = bed[BED6_FIELDS]
    return bed


@click.command()
@click.argument("gtf")
@click.argument("out_bed")
@click.option("--name", "-n",
    default="gene_name",
    show_default=True,
    help="Use which section as bed name, for example 'gene_name', 'gene_id', 'gene_version' ...")
def gtf2bed(gtf, out_bed, name):
    df_gtf = read_gtf(gtf)
    df_bed = gtfdf2beddf(df_gtf, name)
    df_bed.to_csv(out_bed, sep="\t", header=False, index=False)


if __name__ == "__main__":
    eval("gtf2bed()")
