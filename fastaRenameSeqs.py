#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    fastaRenameSeqs.py <fasta> <map> > renamed.fasta

    Description:
    ------------
    Renames sequence headers in a FASTA file per the two-column 
    tab-delimited mapping file <map> where the first column contains 
    the old name and the second column contains the new name.
        ''')
    sys.exit(0)


# output help information if missing command line arguments
args = sys.argv
if len(args) < 3:
    help()
# read mapping file into a dictionary 
map_dct = {}
with open(args[2]) as in_fl:
    for line in in_fl:
        old, new = line.strip().split('\t')
        map_dct[old] = new
# read fasta file and modify each sequence header according to the
# mapping file. if a sequence header is not found in the mapping file
# it is reported on stderr and the original header is output unchanged
with open(args[1]) as in_fl:
    for line in in_fl:
        if line.startswith('>'):
            old = line.strip().lstrip('>')
            try:
                new = '>{0}'.format(map_dct[old])
                print(new)
            except KeyError:
                print('WARNING: Name not found in mapping file\t{0}'.format(
                                                         old), file=sys.stderr)
                print(line, end='')
        else:
            print(line, end='')
