#!/usr/bin/env python3

import sys
from gff3line import GFF3_line


def help():
    print('''
    Usage:
    ------------
    gffSubset -list <path> -gff <path>

    Description: 
    ------------
    Takes a list of LTR_retrotransposon IDs and a GFF3 file and outputs
    feature blocks for the elements in the list. Assumes the GFF3 is
    properly sorted with parent-child relationships preserved.
        ''', file=sys.stderr)
    sys.exit(0)


args = sys.argv
# output help information if missing command line arguments
if ('-list' not in args 
     or '-gff' not in args
     or len(args) < 5):
    help()
# read command line arguments
gff = args[args.index('-gff') + 1]
lst = args[args.index('-list') + 1]
# read list contents into a set
lstSet = set()
with open(lst) as fl:
    for line in fl:
        lstSet.add(line.strip().lstrip('LTR_retrotransposon'))
# read gff lines and output lines belonging to each block whose
# repeat_region feature number is in lstSet
with open(gff) as fl:
    output = None
    for line in fl:
        # skip commented lines
        if line.startswith('#'):
            continue
        GFF = GFF3_line(line)
        if GFF.type == 'repeat_region':
            output = False
            if GFF.attributes['ID'].lstrip('repeat_region') in lstSet:
                output = True
                print('###')
                print(line.strip())
        elif output == True:
            print(line.strip())
