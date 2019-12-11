#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    blastBestHit.py < blast.table > output.table
    
    Description:
    ------------
    This script takes a blast output format 7 table on stdin and outputs 
    the line containing the first hit, which, if the input table is
    unaltered blastn, blastp, etc. output will be the highest-scoring
    hit.
        ''')
    sys.exit(0)

# output help information if requested
if '-h' in sys.argv or '-help' in sys.argv or '--help' in sys.argv:
    help()
# a string that is very unlikely to be the name of a queried sequence
# to serve as the "last query" for the first iteration of the loop below
last_query = '!@#$%^&*()_QWERTYUIOP{'
# for each line in the input blast table output the first hit which if
# the input is unaltered blastn, blastp, etc. output is the hit with the
# highest score
for line in sys.stdin:
	if not line.startswith('#') and not line.startswith('{0}\t'.format(
                                                                  last_query)):
		print(line, end='')
		last_query = line.strip().split()[0]

