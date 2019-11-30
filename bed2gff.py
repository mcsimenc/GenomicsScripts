#!/usr/bin/env python3
# Converts a 10-column BED file to a GFF3 file

import sys


def help():
    print('''
    Usage:
    ------------
    bed2gff < file.bed > file.gff

    Description:
    ------------
    Converts a 10-column BED file such as the type created from a GFF3
    file by the BEDOPS script gff2bed.

    BED uses 0-based numbering for features. Length = start - end
    GFF3 uses 1-based numbering for features. Length = start - end + 1

        ''')
    exit(0)


# Print help information
args = sys.argv
if '-h' in args:
    help()

# Read from stdin
for line in sys.stdin:
    # Output comments from the input BED as they are
    if line.startswith('#'):
        print(line)
    # Read each line in the BED file, convert the coordinates from
    # 0-based to 1-based and output GFF3 lines in which the fields are
    # mostly the same but ordered differently
    else:
        (scaf, bed_start, bed_end, name, score, strand, source, Type, unknown, 
                attr)  = line.strip().split('\t')

        gff_start = str(int(bed_start) + 1)
        gff_end = bed_end

        print('\t'.join([scaf, source, Type, gff_start, gff_end, score, 
                strand, '.', attr]))
