#!/usr/bin/env python3

import sys
import re

def help():
    print('''
    Usage:
    ------------
    gffFilter.py [options] < input.gff > output.gff

    Description:
    ------------
    Removes lines from a GFF3 file if the lines do or do not match
    scaffold names or ID attributes present in lists provided with
    either -scaf or -id.

    Options:
    ------------
    -scaf <path>    Output lines where the scaffold name is present in
                    the list in a file provided

    -id <path>      Output lines where the value of the ID attribute
                    matches an item of the list in the file provided

    -v              Invert. Print non-matching lines.
    ''')
    sys.exit(0)


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


# output help information if not enough arguments were proivded
args = sys.argv
if '-h' in args or '-help' in args or len(args) < 3:
    help()
# parse command line arguments
if '-id' in args:
    IDpth = args[args.index('-id') + 1]
    IDset = set()
    ID = True
    with open(IDpth) as IDfl:
        for line in IDfl:
            IDset.add(line.strip())
else:
    ID = False
if '-scaf' in args:
    Scafpth = args[args.index('-scaf') + 1]
    Scafset = set()
    SCAF = True

    with open(Scafpth) as Scaffl:
        for line in Scaffl:
            Scafset.add(line.strip())
else:
    SCAF = False
if SCAF and ID:
    print(('You used both -scaf and -id but only -id will be used. '
                        'No implementation exists for both.'), file=sys.stderr)
# for each line in the input gff3 file, read the line as a GFF3_line
# object and output lines requested per command line arguments
for line in sys.stdin:
    # output commented lines
    if line.startswith('#'):
        print(line, end='')
        continue
    gffLine = GFF3_line(line)
    if ID:
        if gffLine.attributes['ID'] in IDset:
            if '-v' in args:
                continue
            else:
                print(line, end='')
        else:
            if '-v' in args:
                print(line, end='')
            else:
                continue
    elif SCAF:
        if gffLine.seqid in Scafset:
            if '-v' in args:
                continue
            else:
                print(line, end='')
        else:
            if '-v' in args:
                print(line, end='')
            else:
                continue
