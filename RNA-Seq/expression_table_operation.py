import click
import pandas as pd


GTF_FIELDS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes", "comments"]


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


def read_expression_table(path):
    skip = infer_skip_header_rows(path)
    df = pd.read_table(path, sep="\t", index_col=0, skiprows=skip)
    return df


def gene_id2gene_name(ids, ref_ids, ref_names):
    df = pd.DataFrame({"name": ref_names, "ids": ref_ids})
    df.index = df['ids']
    del df['ids']
    names = df.loc[ids]['name']
    return names


@click.group()
def exp_tab_op():
    pass


@exp_tab_op.command()
@click.argument("table")
@click.argument("gtf")
@click.argument("output")
def index2gene_name(table, gtf, output):
    """
    Change the expression index to gene name.

    \b
    Parameters
    ----------
    table : str
        Path to input tab split table.
    gtf : str
        Path to reference gtf file.
    output : str
        Path to output expression file.
    """
    df_gtf = read_gtf(gtf)
    df_exp = read_expression_table(table)
    gene_names = extract_section_from_attributes("gene_name", df_gtf.attributes)
    gene_ids = extract_section_from_attributes("gene_id", df_gtf.attributes)
    names = gene_id2gene_name(df_exp.index, gene_ids, gene_names)
    df_exp.index = names
    df_exp.to_csv(output, sep="\t")


if __name__ == "__main__":
    eval("exp_tab_op()")
