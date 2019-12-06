#!/usr/bin/env python3

import sys


def help():
    print("""
    Usage:
    ------------
     mergeGFF-Overlaps.py < file.gff > output.gff

    Description:
    ------------
    Combines overlapping features in the input GFF3 file into a single 
    feature. Score and attributes are removed. Commented lines are
    ignored and not output.

    Options:
    ------------
    -keepAttr    For features that do not overlap retain the
                 attributes column from the original GFF3.
""", file=sys.stderr)
    sys.exit()


class GFF3_line:
    """A class to represet GFF3 lines and allow modification of the
    values of its fields.

    Attributes:
    ------------
    field0, ..., field8   strings containing the values of each field
                          in a GFF3 line
    attributes            a dictionary containing the key-value pairs
                          in the GFF3 line 9th field

    Methods:    
    ------------
    str()               Outputs GFF3 line
    repr()              Outputs GFF3 line
    refreshAttrStr()    This needs to be called if changes were made to
                        any of the attributes. It refreshes
    """

    def __init__(self, line):
        """GFF3_line is initialized to contain the fields of the GFF3
        line provided as attributes. The attributes are kept in a
        dictionary and a the keys are ordered in a list to preserve
        the order of attributes upon getting the sring representation
        of the line from the GFF3_line object.
        """
        (self.seqid, 
         self.source, 
         self.type, 
         self.start, 
         self.end, 
         self.score, 
         self.strand, 
         self.phase, 
         self.attributes_str) = line.strip().split('\t')
        # preserve attribute order as a list of keys (attributes_order)
        attributes_list = self.attributes_str.split(';')
        self.attributes_order = [attr.split('=')[0] for attr in 
                                                               attributes_list]
        # store attribute keys and their values in a dictionary
        self.attributes = {attr.split('=')[0]:attr.split('=')[1] for attr in 
                                                               attributes_list}
        # rename the name attribute key to Name so it conforms to the
        # GFF3 specification, where Name is a reserved attribute key
        if 'name' in self.attributes:
            self.attributes['Name'] = self.attributes.pop('name')
            self.attributes_order[self.attributes_order.index('name')] = 'Name'

    def __repr__(self):
        """Output for overloaded functions str() and repr()"""
        return '\t'.join([str(self.seqid), 
                          str(self.source), 
                          str(self.type), 
                          str(self.start), 
                          str(self.end), 
                          str(self.score), 
                          str(self.strand), 
                          str(self.phase), 
                          str(self.attributes_str)])

    def refreshAttrStr(self):
        """If the attributes dictionary or attributes_order has been 
        altered this should be called to update attributes_str.
        """
        self.attributes_str = ';'.join(['='.join(
             [attr, self.attributes[attr]]) for attr in self.attributes_order])


def mergeCoords(A, B):
    """Takes two tuples of coordinates A and B as arguments; A must
    have a start coordinate before or equal to the start coordinate
    of B. If they do not overlap then A and B are returned as input
    and if they overlap the minimum and maximum values are returned.

    let A = (a1, a2), B = (b1, b2) | a1<=b1, a1<=a2, b1<=b2

    case 1: a2<=b1 ---> output A and B

    case 2: b1<a2 && b2>a2 ---> output (a1, b2)

    case 3: b2<=a2 ---> output A
    """

    assert min(A) <= min(B), ("tuples given to mergeCoords in wrong order: "
                                "A={0}, B={1}").format(A,B)
    if min(B) >= max(A):
        return ((A, B), 0)
    elif min(B) < max(A) and max(B) > max(A):
        output = (min(A),max(B))
        return ((output, output), 1)
    elif max(B) <= max(A):
        return ((A, A), 2)
    else:
        raise Exception(("Unexpected result from mergeCoords(A,B) using "
                           " A={0}, B={1}").format(A,B))


# print help information if asked for
if '-h' in sys.argv:
    help()
# for each line in the input gff store as a GFF3_line object in a
# dictionary organized by scaffold. then merge features on each
# scaffold and output
gffFeats = {}
for line in sys.stdin:
    # skip comment lines
    if not line.startswith('#'):
        lineObj = GFF3_line(line)
        scaf = lineObj.seqid
        if scaf in gffFeats:
            gffFeats[scaf].append(lineObj)
        else:
            gffFeats[scaf] = [lineObj]
# for each scaffold sort features by start coordinate then merge 
# overlapping features
for scaf in gffFeats: 
    gffFeats[scaf] = sorted(gffFeats[scaf], key=lambda x:int(x.start))
    newScafFeats = []
    currentFeat = gffFeats[scaf][0]
    i=0
    while i < len(gffFeats[scaf]) - 1:
        nextFeat = gffFeats[scaf][i + 1]
        mergeResult = mergeCoords((int(currentFeat.start), int(currentFeat.end)), 
                                  (int(nextFeat.start), int(nextFeat.end)))
        # feats do not overlap
        if mergeResult[1] == 0:
            currentFeat.start, currentFeat.end = mergeResult[0][0]
            currentFeat.score = '.'
            if '-keepAttr' not in sys.argv:
                currentFeat.attributes_str = '.'
            newScafFeats.append(currentFeat)
            currentFeat = nextFeat
        # feats overlap. continue iterations and check for overlap with
        # subsequent feature
        else:
            currentFeat.start, currentFeat.end = mergeResult[0][0]
        i += 1
    # finish processing last feature
    currentFeat.score = '.'
    if '-keepAttr' not in sys.argv:
        currentFeat.attributes_str = '.'
    newScafFeats.append(currentFeat)
    # replace existing 
    gffFeats[scaf] = newScafFeats
# output new gff lines in order of scaffold and start coordinate
for scaf in sorted(gffFeats.keys()):
    for line in gffFeats[scaf]:
        print(str(line))
