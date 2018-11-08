import sys
import gzip
import subprocess
import logging

import fire

log = logging.getLogger(__name__)
log.addHandler(logging.StreamHandler(stream=sys.stdout))
log.setLevel(level=logging.DEBUG)

def pairix_query(pairs, region1, region2):
    """
    warp of pairix query
    """
    cmd = "pairix {} '{}|{}'".format(pairs, region1, region2)
    #log.debug(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    result, _ = p.communicate()
    #log.debug(result)
    all_items = [tuple(line.split()) for line in result.split("\n")]
    return all_items

def loop_diff(pairs1, pairs2, common, specific1, specific2, distance=10000):
    """
    find the different between loop pairs1 and loop pairs2
    pairs1 and pairs2 must be indexed use pairix.

    :pairs1: (str) path to indexed gz file.
    :pairs2: 
    :common: output file for record common loop
    :specific1: record specific in pairs1
    :specific2: record specific in pairs2
    :distance: two loop consider as same loop,
               if the distance between them less than this.
    """
    common_set = set()
    spec_1 = set()
    spec_2 = set()

    with gzip.open(pairs1) as f, open(specific1, 'w') as fw1:
        for line in f:
            items_1 = tuple(line.strip().split())
            _, chr1, c1, chr2, c2 = items_1[:5]
            c1, c2 = int(c1), int(c2)
            region1 = chr1+':'+str(c1-distance//2)+'-'+str(c1+distance//2)
            region2 = chr2+':'+str(c2-distance//2)+'-'+str(c2+distance//2)
            in_pairs2 = pairix_query(pairs2, region1, region2)
            if in_pairs2 != [()]:
                for items_2 in in_pairs2:
                    common_set.add(items_2)
                    common_set.add(items_1)
            else:
                fw1.write("\t".join(items_1) + "\n")

    with gzip.open(pairs2) as f, open(specific2, 'w') as fw2:
        for line in f:
            items = tuple(line.strip().split())
            if items not in common_set:
                fw2.write("\t".join(items) + "\n")

    with open(common, 'w') as f:
        for items in common_set:
            f.write("\t".join(items) + "\n")


if __name__ == "__main__":
    fire.Fire(loop_diff)
