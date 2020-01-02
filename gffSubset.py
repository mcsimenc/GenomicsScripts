#!/usr/bin/env python3

import sys
from gff3line import GFF3_line


def help():
    print('''
    Usage:
    ------------
    gffSubset -attr <str> -list <path> -gff <path>

    Description: 
    ------------
    Takes a GFF3 file, a list, and the name of an attribute key present
    in the GFF3 file and outputs each line for which the value of the
    attribute key is in the list.
        ''', file=sys.stderr)
    sys.exit(0)


args = sys.argv
# output help information if missing command line arguments
if ('-attr' not in args 
     or '-list' not in args 
     or '-gff' not in args
     or len(args) < 7):
    help()
# read command line arguments
attr = args[args.index('-attr') + 1]
gff = args[args.index('-gff') + 1]
lst = args[args.index('-list') + 1]
# read list contents into a set
lstSet = set()
with open(lst) as fl:
    for line in fl:
        lstSet.add(line.strip())
# read gff lines and output lines containing a key-value pair of the
# attribute given by -attr and one of the elements in the lstSet
with open(gff) as fl:
    for line in fl:
        # skip commented lines
        if line.startswith('#'):
            continue
        GFF = GFF3_line(line)
        if GFF.attributes[attr] in lstSet:
            print(line.strip())
