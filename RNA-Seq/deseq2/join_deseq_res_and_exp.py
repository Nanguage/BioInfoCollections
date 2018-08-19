help = """
Join deseq2 results table with the expression values.

Usage: python join_deseq_res_and_exp.py DESEQ_RES EXPRESS_TABLE COLS

Parameters:
    DESEQ_RES: deseq2 result table
    EXPRESS_TABLE: expression table
    COLS: The column index numbers(start from 0) of expression table,
          for example: 2,3,4 

Example:
    $ python join_deseq_res_and_exp.py ./deseq2_result.tsv ./expression_tab.tsv 1,2,3,4,5,6 > result.tsv

"""

import sys
import pandas as pd

def read_table(path):
    return pd.read_table(path, index_col=0)

def join_table(tab1, tab2):
    result = tab1.join(tab2)
    return result

def select_table(table, cols):
    cols = list(map(int, cols.split(",")))
    return table.iloc[:, cols]

def compose(deseq_tab, exp_tab, cols):
    deseq = read_table(deseq_tab)
    exp = read_table(exp_tab)
    sub_exp = select_table(exp, cols) 
    joined = join_table(deseq, sub_exp)
    result = joined
    return result

def output(result, file=sys.stdout, sep='\t'):
    if not isinstance(file, str):
        result.to_csv(file, sep=sep)
        return
    with open(file, 'w') as buf:
        result.to_csv(buf, sep=sep)

def deal_help(args):
    if len(args) != 4 or\
       '-h' in args or\
       '--help' in args:
        print(help)
        sys.exit(0)    
    
if __name__ == "__main__":
    args = sys.argv
    deal_help(args)
    
    deseq_table = args[1]
    exp_table = args[2]
    cols = args[3]

    result = compose(deseq_table, exp_table, cols)
    output(result)
