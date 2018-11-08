import os
import subprocess

import fire

def looplist2pairs(in_file, out_file, index=True):
    """
    convert loop list to pairs format and index it.
    """
    with open(in_file) as fi, open(out_file, 'w') as fo:
        for i, line in enumerate(fi):
            if i == 0:
                continue
            items = line.strip().split()
            chr1, start1, end1, chr2, start2, end2 = items[:6]
            name = os.path.basename(in_file) + '_' + str(i)
            start1, end1, start2, end2 = [int(i) for i in start1, end1, start2, end2]
            c1, c2 = (start1 + end1)//2, (start2 + end2)//2
            # to upper trangle
            if chr1 > chr2:
                chr1, c1, chr2, c2 = chr2, c2, chr1, c1
            if chr1 == chr2 and c1 > c2:
                c1, c2 = c2, c1
            line = '\t'.join([name, chr1, str(c1), chr2, str(c2), '+', '+'])
            fo.write(line + "\n")

    # sort pairs file
    subprocess.check_call(
        "sort -k2,2 -k3,3n -k4,4 -k5,5n {} > {}".format(out_file, out_file+'.sorted'),
        shell=True)

    # index
    subprocess.check_call(
        ['mv', out_file+'.sorted', out_file]
    )

    # index
    subprocess.check_call(['bgzip', '-f', out_file])
    subprocess.check_call(['pairix', out_file+'.gz'])


if __name__ == "__main__":
    fire.Fire(looplist2pairs)
