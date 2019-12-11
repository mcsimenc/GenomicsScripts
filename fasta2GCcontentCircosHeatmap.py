#!/usr/bin/env python3

# If you have questions about this script or something doesn't work right you can email Matt at mcsimenc@gmail.com

import sys
from math import ceil

def help():
    print('''
    Usage:
    ------------
    fasta2GCcontentHeatmapTrack.py -window <int> \\
                                   -fasta <path> \\
                                   [-scafList <path>]

    Description:
    ------------
    From a FASTA file, GC content per -window bp is calculated and
    output in Circos heatmap track format. Optionally restrict output
    to a list of sequences with -scafList.

    Options:
    ------------
    -fasta      <path>    FASTA file
    -window     <int>     Bin size in bp
    -scafList   <path>    A list of sequences to limit output to

    Output:
    ------------
    4-column tab-delimited file in Circos heatmap track format

        scaf    start    stop    density
    ''', file=sys.stderr)
    sys.exit()


def genWindowCoords(n, windowLen):
    """Generate coordinates for the nth window"""
    return (n*windowLen, (n+1)*windowLen-1)


def seq2gc(seq):
    """Calculates GC content for a sequence"""
    seq = seq.lower()
    return (seq.count('g') + seq.count('c')) / len(seq)


def seq2gcHeatmap(seqName, seq, windowLen):
    """
    From a nucleotide sequence, its name, and a window length, output 
    Circos heatmap track format lines for GC content in the sequence
    for each window.
    """
    seqLen = len(seq)
    numWindows = ceil(seqLen // windowLen)
    # for each window generate its start and end coordinate, calculate
    # gc content, and output a circos heatmap track line
    for i in range(numWindows):
        # generate window coordinates
        window = genWindowCoords(i, windowLen)
        # calculate gc content for this window
        gc = seq2gc(seq[window[0]:window[1]+1])
        print('{0}\t{1}\t{2}\t{3}'.format(seqName, window[0], window[1], gc))
    # calculate and output gc content for the last, potentially smaller
    # window
    window = genWindowCoords(numWindows, windowLen)
    gc = seq2gc(seq[window[0]:seqLen])
    print('{0}\t{1}\t{2}\t{3}'.format(seqName, window[0], seqLen, gc))


# output help information if asked for or command line arguments are
# missing
args = sys.argv
if ('-h' in args 
    or len(args) < 5 
    or '-fasta' not in args 
    or '-window' not in args):
    help()
# read command line arguments
if '-scafList' in args:
    scafListFl = args[args.index('-scafList') +1]
    scafList = open(scafListFl).read().strip().split('\n')
else:
    scafList = None
windowLen = int(args[args.index('-window') +1])
fasta_filepath = args[args.index('-fasta') +1]
# for each sequence in the fasta file, read the sequence into memory
# and output circos heatmap track lines
seq = ''
seqName = ''
with open(fasta_filepath) as fastaFl:
    for line in fastaFl:
        if line.startswith('>'):
            if seq != '':
                # restrict output to sequences in -scafList, if provided
                if scafList != None:
                    if seqName in scafList:
                        seq2gcHeatmap(seqName, seq, windowLen)
                # if -scafList not provided, output for all sequences
                else:
                    seq2gcHeatmap(seqName, seq, windowLen)
            seqName = line.strip()[1:]
            seq = ''
        else:
            seq += line.strip()
    # output for the last sequence
    if scafList != None:
        if seqName in scafList:
            seq2gcHeatmap(seqName, seq, windowLen)
    else:
        seq2gcHeatmap(seqName, seq, windowLen)
