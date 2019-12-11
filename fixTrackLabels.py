#!/usr/bin/env python3

import sys
import os


def help():
    print('''

    Usage:
    ------------
    fixTrackLabels.py -track <file> -karyotype <file>

    Description:
    ------------
    Replaces track labels with integer IDs found in karyotype file.

    Example:
    ------------
    karyotype file:
        chr - 0    scaffold00001_v1.0    0    1291680    0

    track file:
        scaffold00001_v1.0    0    99999    63.2241

    new track file:
        0    0    99999    63.2241
    ''', file=sys.stderr)
    sys.exit()


# output help information if missing command line arguments
args = sys.argv
if not '-track' in args and not '-karyotype' in args:
    help()
# read command line arguments
trackFl = args[args.index('-track') +1]
karyotypeFl = args[args.index('-karyotype') +1]
# read karyotype file and create map of old identifier to new identifier
_map = {}
with open(karyotypeFl, 'r') as inFl:
    for line in inFl:
        _line = line.strip().split()
        _map[_line[3]] = _line[2]
# for each line in the track file, modify accordingly and output new
# line
with open(trackFl, 'r') as inFl:
    for line in inFl:
        _line = line.strip().split()
        _line[0] = _map[_line[0]]
        print('\t'.join(_line))
        
