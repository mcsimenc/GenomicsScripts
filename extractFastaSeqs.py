#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    extract_fasta_seqs.py [fasta_file] [out_file] [options]

    Description:
    ------------
    Outputs a file [out_file] in FASTA format containing the 
    sequences derived from [fasta_file] whose sequence headers
    correspond to the names given by -seq_list.

    Options:
    ------------
    -seq_list    List of sequence identifiers
    -v           Invert: select all except sequences in -seq_list
        ''')
    sys.exit(0)


args = sys.argv
# output help information if not enough command line arguments are 
# provided or if -seq_list is missing
if ("-help" in args or "-h" in args 
  or len(args) < 5 
  or '-seq_list' not in args):
    help()
# read command line arguments
fasta_file = open(args[1])
out_file = open(args[2],"w")
seq_list = open(args[args.index('-seq_list')+1])
# add the sequence names given by -seq_list to a set
seq_set = set()
for seq in seq_list:
    seq_set.add(seq.strip())
# for each line in the input fasta file, mark a sequence as a match
# if the sequence header is in the -seq_list and output or not depending
# on whether -v is given
match = 0
for line in fasta_file:
    if line.startswith(">"):
        match = 0
        seqname = line.strip()[1:]
        if '-v' in args:
            if seqname not in seq_set:
                out_file.write(line)
                match = 1
            else:
                seq_set.remove(seqname)
        elif '-v' not in args:
            if seqname in seq_set:
                match = 1
                out_file.write(line)
                seq_set.remove(seqname)
    else:
        if match == 1:
            out_file.write(line)
        elif match == 0:
            continue
if not seq_set == set():
    print('These seqences were not found:', file=sys.stderr)
    for seq in seq_set:
        print(seq)
