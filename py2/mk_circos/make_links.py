# -*- coding: utf-8 -*-

import os
import re

from stringxmlparse import get_iid2gid, get_eid2type, get_connection
from gtfparser import GtfParser
from model import Link, Label


STRING_XML_FILE = "../data/up.xml"
GTF_FILE = "../../genome/Mycobacterium_tuberculosis_h37ra.ASM1614v1.32.gtf"
LINKS_DIR = "../draw/data/links_up/"
LABEL_FILE = "../draw/data/band_labels_up.txt"


def get_gid2pos(gff_file):
    """Get gene_id 2 gene position information from gtf file."""
    gid2pos = {}
    parser = GtfParser(gff_file)
    for item in parser:
        gene_id  = item.gene_id
        pos      = item.position
        gid2pos[gene_id] = str(pos[0]), str(pos[1])
    return gid2pos


def to_links_files(string_xml_file, out_dir, gtf_file, chr_):
    """Write links to file."""
    connections = get_connection(string_xml_file)
    gid2pos = get_gid2pos(gtf_file)
    type2links = {}
    for gid1, gid2, link_type in connections:
        type2links.setdefault(link_type, [])
        start1, end1 = gid2pos[gid1]
        start2, end2 = gid2pos[gid2]
        link = Link(chr_, start1, end1, chr_, start2, end2, link_type)
        type2links[link_type].append(link)
    for link_type in type2links:
        with open(os.path.join(out_dir, 
            '_'.join(re.split('[/ ]', link_type))+'.txt'), 'w') as f:
            for l in type2links[link_type]:
                f.write(str(l) + '\n')


def to_label_file(string_xml_file, out_file, gtf_file, chr_):
    """Write label file."""
    connections = get_connection(string_xml_file)
    gid2pos = get_gid2pos(gtf_file)
    gids = set()
    for gid1, gid2, link_type in connections:
        gids.add(gid1)
        gids.add(gid2)
    with open(out_file, 'w') as f:
        for gid in gids:
            s, e = gid2pos[gid]
            l = Label(chr_, str(s), str(e), gid)
            f.write(str(l) + '\n')
    

if __name__ == '__main__':
    to_links_files(STRING_XML_FILE, LINKS_DIR, GTF_FILE, 'mtb')
    to_label_file(STRING_XML_FILE, LABEL_FILE, GTF_FILE, 'mtb')
