#!/usr/bin/env python3

from gff3line import GFF3_line
import sys


# Using LTRdigest-output GFF3 as input, calculate the distance of each
# polypurine tract from their nearest long terminal repeat. It is assumed that
# the nearest LTR is the 3'-LTR

# initialize variables before loop
len1, len2, ltr1end, ltr2start, pptstart, pptend = [None] * 6
for line in sys.stdin:
    # skip commented lines
    if line.startswith('#'):
        continue
    # read GFF line into a GFF3_line() type
    gffline = GFF3_line(line)
    # first or second ltr
    if gffline.type == 'long_terminal_repeat':
        # if none, this is the first LTR
        if ltr1end == None:
            ltr1end = int(gffline.end)
        # otherwise this is the second LTR
        else:
            ltr2start = int(gffline.start)
        # pptstart will be None at all first-LTR features and, if present, will
        # not be None at every second LTR
        if pptstart != None:
            len1 = pptstart - ltr1end + 1
            len2 = ltr2start - pptend + 1
            if len1 < len2:
                print(len1)
            else:
                print(len2)
        # whether or not RR_tract was encountered in this element, reinitialize
        # variables to None because this marks the end of an element
        if ltr2start != None:
            len1, len2, ltr1end, ltr2start, pptstart, pptend = [None] * 6
            
    # get both start and end for the PPT
    elif gffline.type == 'RR_tract':
        pptstart = int(gffline.start)
        pptend = int(gffline.end)
