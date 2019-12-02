#!/usr/bin/env python3

# If you have questions about this script or something doesn't work right you can email Matt at mcsimenc@gmail.com

import sys


def help():
    print('''

    Usage:
    ------------
    gff2circos-tile.py -gff <file> [-valueDef <file>]

    Description:
    ------------
    Takes a GFF3 file as input and outputs a tile track file for
    Circos. The fourth column in the output (value) tells Circos 
    how to color the feature just like it does with heatmap tracks. 
    The default value is 0 but a two-column file can be provided using
    valueDef to customize values.

    Required parameters:
    ------------
    -gff          <path>    Path to input maker gff3 file (mandatory).

    -valueDef     <path>    Tab-delimited file where first column is a 
                            string to search for in each GFF3 line and 
                            the second column is the value to assign if 
                            the string is found. The second column needs 
                            to contain numbers only. e.g: exon  1
                                                          gene  2
                                                          intron    3

    Output:
    ------------
    scaf    start    stop    value
''', file=sys.stderr)


def gff2circosTileTrack(gffFl, valueDef=None):
    """Reads a GFF3 file and outputs a line in Circos tile-track format
    for each feature.
    """
    # read gff file
    with open(gffFl) as fl:
        for line in fl:
            # skip commented lines
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                seq = fields[0]
                start = fields[3]
                end = fields[4]
                value = 0
                # read the value definition file if provided, search
                # each gff line for a match to one of the definitions.
                # if multiple definition file lines match a single gff
                # line, report it to stderr
                if not valueDef == None:
                    foundDef = False
                    MatchValue = 0
                    for definition in valueDef:
                        if definition[0] in line:
                            value = definition[1]
                            if foundDef == True and value != MatchValue:
                                print(('Found more than one definition with '
                                'different values for GFF3 line:\n{0}').format(
                                                line.strip()), file=sys.stderr)
                            foundDef = True
                            MatchValue = value
                print('{0}\t{1}\t{2}\t{3}'.format(seq, start, end, value))

                
if __name__ == '__main__':
    args = sys.argv
    # print help information if not enough or the right gff command line
    # arguments are not provided
    if '-h' in args or len(args) < 3 or '-gff' not in args:
        help()
        sys.exit()
    # save gff file path
    gff_filepath = args[args.index('-gff') +1]
    # read value definitions file
    if '-valueDef' in args:
        valueDefFl = args[args.index('-valueDef') +1]
        valueDef = set([(definition.split('\t')[0], 
                         float(definition.split('\t')[1])) for definition in 
                                 open(valueDefFl).read().strip().split('\n') ])
    else:
        valueDef = None
    # convert gff to Circos tile track format
    gff2circosTileTrack(gff_filepath, valueDef)
