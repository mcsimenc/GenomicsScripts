#!/usr/bin/env python3

from Bio import SeqIO
import sys


def help():
    print('''
    Usage:
    ------------
    fasta2circosKaryotype.py -f <path> [-c]

    Description:
    ------------
    Outputs the lengths of the sequences in the input FASTA file in
    Circos karyotype format for drawing ideograms. Sequences are ordered
    descending by length.

    Options:
    ------------
    -c      Use if sequence names are Chr01, Chr02, ...
    ''')
    sys.exit()


def fasta2track(fasta):
    """
    Returns a list of tuples (seqname, len(seq)) for seqs in fasta
    """
    # read sequences and sort by length
    seq_lens = sorted([(seq.name, len(seq)) for seq 
                    in SeqIO.parse(fasta, format='fasta')], key=lambda x:x[0])
    # unless sequences have implicit ordering based on their names,
    # print in descending order by length
    for i, seq in enumerate(seq_lens):
        if '-c' in args:
            chrNum = int(seq[0].lstrip('Chr'))
        else:
            chrNum = i
        print('chr - {0}\t{1}\t0\t{2}\t{3}'.format(chrNum, seq[0], seq[1], i))


if __name__ == '__main__':
    # output help information if missing command line arguments
    args = sys.argv
    if '-h' in args or '-f' not in args:
        help()
    # read command line arguments
    fastaPth = args[args.index('-f')+1]
    # convert to circos track format and output
    fasta2track(fastaPth)
