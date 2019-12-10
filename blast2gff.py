#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    blast2gff.py -blast <file> > <output.gff>

    Description:
    ------------
    Converts blastn, blastp, etc. tabular output (-outfmt 6 or 7)
    to GFF3 format.

    Blast formats 6 and 7 can be customized so as to include
    subject strand, which is assumed to be the 13th field in the 
    blast file (if included). If subject strand is not in the blast
    output then strand is . (period) in the output GFF3. 

    Options:
    ------------
    -map    identifier map from which to obtain a value for the GFF3
            ID attribute  e.g. for NCBI nt or nr. Before first space
            is expected to be the identifier and after the space
            is a description. For example:
                 X17276.1 Giant Panda satellite 1 DNA
        ''')
    sys.exit(0)


# output help information if command line arguments are missing
args = sys.argv
if not '-blast' in args or len(args) < 3:
    help()
# read command line arguments
blast_fl_pth = args[args.index('-blast') + 1]
idMap = {}
if '-map' in args:
    with open(args[args.index('-map') + 1]) as mapFl:
        for line in mapFl:
            item = line.strip().split(' ', 1)
            idMap[item[0]]=item[1]
# output the top line in a GFF3 file
print('##gff-version 3')
# these fields are the same for every line in the output gff
source = "blastn"
type = "hit"
phase = "."
# for each line in the input blast table
with open(blast_fl_pth) as blast_fl:
    for line in blast_fl:
        # ignore commented lines such as those in blast -outfmt 7
        if line.startswith("#"):
            continue
        else:
            blast_line = line.strip().split()
            # if -map is used, obtain the attributes field ID value from
            # the -map, otherwise use the identifier in the first field
            # of the blast table as the value for the ID attribute
            hit = blast_line[1]
            if '-map' in args:
                try:
                    desc = idMap[hit]
                    if ' ' in desc:
                        desc = desc.replace(' ', '_')
                    attr = "ID={0}_{1}".format(hit, desc)
                except KeyError:
                    print('Not found in map file\t{0}'.format(hit), 
                                                               file=sys.stderr)
            else:
                attr = "ID={0}".format(hit)
            # read sequence name, start, end, and score
            seq_id = blast_line[0]
            start = int(blast_line[6])
            end = int(blast_line[7])
            score = blast_line[10]
            # if hit is on reverse strand swap start and end
            if start > end:
                start, end = end, start
            # attempt to read strand information
            try:
                if blast_line[12] == "plus":
                    strand = "+"
                elif blast_line[12] == "minus":
                    strand = "-"
            except IndexError:
                strand = "."
            # output gff lines
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(seq_id, 
                                                                       source, 
                                                                       type, 
                                                                       start, 
                                                                       end, 
                                                                       score, 
                                                                       strand, 
                                                                       phase, 
                                                                       attr))
