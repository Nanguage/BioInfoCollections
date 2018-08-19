import os
import sys
import argparse as argp

fold_change_cutoff = 0.58
pvalue_cutoff = 0.05

LOG2FC_COL_INDEX = 2
PVALUE_COL_INDEX = 5
PADJ_COL_INDEX = 6
BASEMEAN_COL_INDEX = 1

def argument_parser():
    parser = argp.ArgumentParser(
        description="Filter the deffiential analyze result table.")
    parser.add_argument("input", help="input table.")
    parser.add_argument("--pvalue-cutoff", "-p",
        default=0.05,
        type=float,
        help="The Pvalue cutoff, record will be filter out if pvalue large than this.")
    parser.add_argument("--log2-fold-change-cutoff", "-f",
        default=1,
        type=float,
        help="The log2(fold_change) cutoff, record will be filter out if abs(log2fc) small than this")
    parser.add_argument("--mean-cutoff", "-m",
        type=float,
        help="The base mean cutoff, record will be filter out if baseMean small than this")
    parser.add_argument("--padj",
        action="store_true",
        help="Use adjusted pvalue as pvalue.")
    return parser

def parse_argument():
    parser = argument_parser()
    args = parser.parse_args()
    return args


def filter_file(args):
    with open(args.input) as f:
        for i, line in enumerate(f):
            line = line.strip()
            items = line.split()
            if i == 0:
                print(line)
                continue
            log2fc = float(items[LOG2FC_COL_INDEX])

            try:
                idx = PVALUE_COL_INDEX if args.padj is None else PADJ_COL_INDEX
                p = float(items[idx])
            except ValueError:
                continue

            condition = (abs(log2fc) >= args.log2_fold_change_cutoff) and (p <= args.pvalue_cutoff)
            if args.mean_cutoff:
                base_mean = float(items[BASEMEAN_COL_INDEX])
                condition = condition and (base_mean >= args.mean_cutoff)
            if condition:
                yield line + "\n"

def output(result, file=sys.stdout):
    if not isinstance(file, str):
        for line in result:
            file.write(line)
    else:
        with open(file, 'w') as buf:
            for line in result:
                buf.write(line)


if __name__ == "__main__":
    args = parse_argument()
    output(filter_file(args))
