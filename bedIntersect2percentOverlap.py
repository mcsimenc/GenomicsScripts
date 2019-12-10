#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    bedIntersect2percentOverlap.py < <bed_file>

    Description:
    ------------
    Calculates the percent overlap for each feature in a BED format
    file created using intersectBed with the -wa and -wb flags.

    Seven column tab-delimited output:

    1. bp_overlap 
    2. feat1 
    3. feat1Len 
    4. feat1propOverlap 
    5. feat2 
    6. feat2Len 
    7. feat2propOverlap
    ''')
    sys.exit(0)


# output help information if asked for
args = sys.argv
if '-h' in args:
    help()
# output header line for output table
print('bp_overlap\tfeat1\tfeat1Len\tfeat1propOverlap\tfeat2\tfeat2Len\tfeat2propOverlap')
for line in sys.stdin:
    # read bed input
    contents = line.strip().split('\t')
    feat1start = float(contents[1])
    feat1end = float(contents[2])
    feat1name = contents[3]
    feat2start = float(contents[11])
    feat2end = float(contents[12])
    feat2name = contents[13]
    # calculate feature 1 and 2 lengths
    feat1length = abs(feat1end - feat1start)
    feat2length = abs(feat2end - feat2start)
    # calculate bp overlap
    coords = [feat1start, feat1end, feat2start, feat2end]
    coords.pop(coords.index(max(coords)))
    coords.pop(coords.index(min(coords)))
    bp_overlap = abs(coords[0] - coords[1])
    # calculate percent overlap
    percent_overlap1 = bp_overlap/feat1length
    percent_overlap2 = bp_overlap/feat2length
    # output table with the longer feature first
    lengths = [feat1length, feat2length]
    longer = lengths.index(max(lengths))
    if longer == 0:
        print('{0}\t{1}\t{2}\t{3:.3f}\t{4}\t{5}\t{6:.3f}'.format(
                                                             int(bp_overlap), 
                                                             feat1name, 
                                                             int(feat1length), 
                                                             percent_overlap1, 
                                                             feat2name, 
                                                             int(feat2length), 
                                                             percent_overlap2))
    elif longer == 1:
        print('{0}\t{1}\t{2}\t{3:.3f}\t{4}\t{5}\t{6:.3f}'.format(
                                                             int(bp_overlap), 
                                                             feat2name, 
                                                             int(feat2length), 
                                                             percent_overlap2, 
                                                             feat1name, 
                                                             int(feat1length), 
                                                             percent_overlap1))
