"""
Tools for process matrix file form the deeptools `computeMatrix` command.
"""

import gzip

import numpy as np
import click
from deeptools.heatmapper import heatmapper


@click.group()
def deepmatrix_process():
    pass


@deepmatrix_process.command()
@click.argument("input")
@click.argument("output")
def mean(input, output):
    """
    Get the mean value of the matrix,
    for reduce it's size speed up the loading process
    when plot the profile.

    \b
    Args
    ----
    input : str
        Path to input gziped matrix.
    output : str
        Path to output gziped matrix(vector).
    """
    hm = heatmapper()
    hm.read_matrix_file(input)
    mat = hm.matrix.matrix
    mean_mat = np.array(np.nanmean(mat, axis=0)).reshape(1, -1)
    outfile = gzip.open(output, 'wt')
    metainfo = str(hm.parameters)
    outfile.write("@"+metainfo+"\n")
    np.savetxt(outfile, mean_mat)
    outfile.close()


if __name__ == "__main__":
    deepmatrix_process()
