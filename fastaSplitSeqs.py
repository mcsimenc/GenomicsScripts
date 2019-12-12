#!/usr/bin/env python3

import sys


def help():
	print('''
		Usage:
		------------
		fastaSplitSeqs.py <fasta_file> <output_prefix>

		Description:
		------------
        Writes a separate FASTA file for each sequence in the input
        FASTA.
		''')
	sys.exit(0)


args = sys.argv
# output help information if missing command line arguments
if '-help' in args or '-h' in args or len(args) < 2:
	help()
# read command line arguments
in_fl = open(args[1])
try:
	out_prefix = args[2]
except IndexError:
	out_prefix = in_fl
seq = 0
# write a new fasta file for each sequence in the original fasta file
for line in in_fl:
	if line.startswith('>'):
		seq += 1
		out_fl = open(out_prefix + "_" + str(seq) + ".fasta", "w")
		out_fl.write(line)
	else:
		out_fl.write(line)
out_fl.close()
