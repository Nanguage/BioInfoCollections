import click
from pyfaidx import Fasta


def parse_genome_region(region):
    import re
    chr, s, e = re.split("[:-]", region)
    s = int(s)
    e = int(e)
    assert e >= s, "genome region end position should greater than start position."
    return (chr, s, e)


def probe_regions(chrom, start, end, probe_len, probe_num):
    region_len = end - start
    block_len = region_len // probe_num
    s_ = start
    for _ in range(probe_num):
        e_ = s_ + probe_len
        r = (chrom, s_, e_)
        yield r
        s_ += block_len

def extract_seq(region_iter, fasta):
    for (c, s, e) in region_iter:
        seq = fasta[c][s:e]
        res = ((c, s, e), seq)
        yield res

def output(res_iter, outpath):
    with open(outpath, 'w') as f:
        for (c, s, e), seq in res_iter:
            f.write(">{}:{}-{}".format(c, s, e) + "\n")
            f.write(str(seq)+"\n")


@click.command()
@click.argument("fasta-path")
@click.argument("output-path")
@click.argument("genome-region")
@click.option("--probe-number", "-n", type=int, default=100,
    help="How many probe. default 100")
@click.option("--probe-length", "-l", type=int, default=30,
    help="How long is one probe. default 30")
def main(fasta_path, output_path, genome_region, probe_number, probe_length):
    """
    \b
    Args
    ----
    fasta_path : str
        Path to fasta file,
        fasta should be indexed with `samtools faidx`.
    output_path : str
        Path to output fasta file.
    genome_region : str
        For example, chr9:8900000-9100000
    """
    print("fasta file: " + fasta_path)
    print("genome region: " + genome_region)
    print("probe number: " + str(probe_number))
    fa = Fasta(fasta_path)
    chrom, start, end = parse_genome_region(genome_region)
    region_length = end - start
    if probe_length * probe_number >= region_length:
        raise ValueError("probe_length({}) * probe_number({}) should less than region length({})".format(probe_length, probe_number, region_length))
    regions = probe_regions(chrom, start, end, probe_length, probe_number)
    results = extract_seq(regions, fa)
    output(results, output_path)


if __name__ == "__main__":
    eval("main()")