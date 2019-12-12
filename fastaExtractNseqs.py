#!/usr/bin/env python3

import sys


def help():
	print('''
		Usage:
		------------
		fastaExtractNseqs.py -fasta <file> -n <int> -j <int>

		Description:
		------------
		Prints the jth n sequences from the FASTA file. -n 2 -j 2 would 
        output the second two sequences.

		Options:
		------------
		-fasta <file>	The input fasta file
		-n <int>	The number of sequences to output
		-j <int>	The set of sequences to output. e.g:
				    -j 1 means output the first n seqs,
				    -j 2 menas output the second n seqs.
		''')
	sys.exit(0)


# output help information if missing command line arguments
args = sys.argv
if (not len(args) == 7 
    or not '-j' in args 
    or not '-n' in args 
    or not '-fasta' in args):
	help()
# read command line argument
in_fl = open(args[args.index('-fasta') + 1])
n = int(args[args.index('-n') + 1])
j = int(args[args.index('-j') + 1])
ct = 0
j_ct = 1
for line in in_fl:
    # increment n at seq headers
	if line.startswith('>'):
		ct += 1
        # increment j if count reaches n
		if ct > n: 
			ct = 1
			j_ct += 1
		if j_ct < j:
			continue
		elif j_ct == j:
			print(line.strip())
        # stop if past the jth set of n
		elif j_ct > j:
			break
	else:
		if j_ct == j:
			print(line.strip())
		else:
			continue
