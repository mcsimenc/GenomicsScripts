#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    gffRemoveScafPart.py [options] < input.gff > output.gff

    Description:
    ------------
    Removes features present in the range provided from a GFF3 file.

    Options:
    ------------
    -scaf  <string>     Scaffold name from which to remove features

    -range <int-int>    Range to remove, e.g. -range 452-1823

    -h                  Output help information
    ''')
    sys.exit(0)


# print help information if not enough command line arguments are
# provided
args = sys.argv
if ('-h' in args 
 or '-help' in args 
 or len(args) < 5
 or '-scaf' not in args
 or '-range' not in args):
    help()
# parse command line arguments
scafTarget = args[args.index('-scaf')+1]
r = args[args.index('-range')+1].split('-')
startTarget = int(r[0])
endTarget = int(r[1])
# for each line in the GFF3 file extract scaffold name and coordinates
# and check if any part of the feature overlaps the region specified
# by -range. If it does, do not output that line. If not, output that
# line.
for line in sys.stdin:
    if not line.startswith('#'):
        contents = line.strip().split('\t')
        scaf = contents[0]
        if scaf == scafTarget:
            start = int(contents[3])
            end = int(contents[4])
            # feature is on scaffold and within range
            if ((start >= startTarget and start < endTarget) 
             or (end > startTarget and end < endTarget)):
                continue
            # feature is on scaffold but out of range
            else:
                print(line, end='')
        # feature is not on scaffold
        else:
            print(line, end='')
