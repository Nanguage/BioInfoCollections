# -*- coding: utf-8 -*-

import re

from gtfparser import GtfParser
from model import Chr, Band


GTF_FILE = "../../genome/Mycobacterium_tuberculosis_h37ra.ASM1614v1.32.gtf"
OUT_FILE = "../draw/data/karyotype.mtb_h37ra.txt"


def get_bands(gtf_file, chr_, colors):
    """ get all bands from gtf file"""
    def rotate_color(colors):
        while True:
            for c in colors:
                yield c
    colors = rotate_color(colors)
    bands = []
    parser = GtfParser(gtf_file)
    for item in parser:
        if item.type != 'exon': continue
        band = Band(chr_.id_, item.gene_id, item.gene_id,
                str(item.position[0]), str(item.position[1]), next(colors))
        bands.append(band)
    return bands


def main():
    chr1 = Chr('mtb', 'mtb', '0', '4419977', 'white')
    colors = ['b2', 'b1', 'b0', 'm', 'r0', 'r1', 'r2']
    colors += reversed(colors[1:-1])
    bands = get_bands(GTF_FILE, chr1, colors)
    with open(OUT_FILE, 'w') as f:
        f.write(str(chr1)+"\n")
        for band in bands:
            f.write(str(band) + "\n")
        

if __name__ == "__main__":
    main()
