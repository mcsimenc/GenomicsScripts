#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    gffRenameScaffolds.py -map <file> -input <file>

    Description:
    ------------
    Rename scaffolds/sequence ids in the GFF3 input file. The file
    given by -map should be a two-column text file where the first
    column contains the old scaffold name and the second column
    contains the new scaffold name.

    Output:
    ------------
    Modified lines from the input GFF3
    ''')
    sys.exit()


# print help information if not enough command line arguments were
# provided
args = sys.argv
if ('-h' in args 
    or len(args) < 5 
    or '-map' not in args 
    or '-input' not in args):
    help()
# read command line arguments
map_flname = args[args.index('-map') + 1]
in_flname = args[args.index('-input') + 1]
map_dct = {}
# for each line put the values of the two columns into a dictionary
# {oldScafName:newScafName}
with open(map_flname) as map_fl:
    for line in map_fl:
        contents = line.strip().split('\t')
        try:
            from_id = contents[0]
            to_id = contents[1]
        except KeyError:
            continue
        if from_id in map_dct:
            print('Duplicate ID in first column: {0}'.format(from_id))
        else:
            map_dct[from_id] = to_id
# for each line in the input gff replace the old scaffold name with the
# new name
with open(in_flname) as input_file:
    for line in input_file:
        # these lines are not gff features but contain scaffold names
        # rename them accordingly
        if  line.startswith('##sequence-region'):
            contents = line.strip().split()
            contents[1] = map_dct[contents[1]]
            print('{0}\t{1}'.format(contents[0], ' '.join(contents[1:])))
        # output commented lines unaltered
        elif line.startswith('#'):
            print(line)
        # replace the scaffold name in the first field of gff feature
        # lines with the new name
        else:
            contents = line.strip().split()
            contents[0] = map_dct[contents[0]]
            print('\t'.join(contents))
