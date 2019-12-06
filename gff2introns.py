#!/usr/bin/env python3

import sys


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


def help():
    print('''
        Usage:
        ------------
        gff2introns.py < input.gff > output.txt

        Description:
        ------------
        Takes a GFF3 on stdin and outputs a list of intron lengths to
        stdout and a list of number of exons for each gene to stderr.
        Instead of a list of lengths and exon counts, GFF3 lines for 
        each intron are output if -gff is used. This script assumes the 
        input GFF3 is sorted by start coordinate and that features of 
        type exon follow the gene of which they are a part.

        Options:
        ------------
        -gff    Output GFF3 of intron features instead of 
        ''')
    sys.exit(0)


# Get parameters
args = sys.argv
# If asked for help or insufficient parameters
if "-help" in args or "-h" in args:
    help()
exon_count = None
last_exon_end = None
gene = None
# read each line of gff input and output one gff3 intron line for each
# pair of successive exons in each gene
if '-gff' in args:
    for line in sys.stdin:
        # ignore commented lines
        if not line.startswith('#'):
            gffLine = GFF3_line(line)
            contents = line.strip().split('\t')
            if contents[2] == 'gene':
                # use exon counts to calculate intron counts for
                # inclusion in the value of the new intron feature gff
                # line's ID attribute
                exon_count = 0
                # reset this variable to mark the first exon in a
                # gene
                last_exon_end = None
                continue
            elif contents[2] == 'exon':
                exon_count += 1
                # this is the second of a pair of successive exons.
                # use the exon coordinates to formulate an intron line
                # with appropriate coordinates
                if not last_exon_end == None:
                    intronStart = last_exon_end + 1
                    intronEnd = int(contents[3]) - 1
                    gffLine.type = 'intron'
                    gffLine.start = str(intronStart)
                    gffLine.end = str(intronEnd)
                    gffLine.score = '.'
                    gffLine.phase = '.'
                    # get the gene name from this exon's ID attribute
                    # and format the intron's ID attribute as
                    # geneName-intron_number
                    # assumes the naming convention used is
                    # geneName-exonName
                    gffLine.attributes['ID'] = gffLine.attributes['ID'].split(
                                 '-')[0] + ':intron_{0}'.format(exon_count - 1)
                    gffLine.refreshAttrStr()
                    # output the new intron gff line
                    print(str(gffLine))
                    last_exon_end = int(contents[4])
                # mark this exon as the first in the gene and continue
                else:
                    last_exon_end = int(contents[4])
# read each line of gff input and output a list of intron lengths and a
# list of exon counts to stdout and sterr respectively
else:
    for line in sys.stdin:
        if not line.startswith('#'):
            contents = line.strip().split('\t')
            if contents[2] == 'gene':
                # output the last gene's exon count to stderr
                if not exon_count == None:
                    print('{0}\t{1}'.format(gene, exon_count), file=sys.stderr)
                exon_count = 0
                last_exon_end = None
                # extract gene ID. assumes the ID key is the first key
                # in the attributes column
                gene = contents[8].split(';')[0].split('=')[1]
                continue
            elif contents[2] == 'exon':
                exon_count += 1
                # if this is the second of a pair of successive exons
                # then calculate the length between them as the intron
                # length and output that to stdout
                if not last_exon_end == None:
                    intronStart = last_exon_end + 1
                    intronEnd = int(contents[3]) - 1
                    intronLen = intronEnd - intronStart + 1
                    print('{0}\t{1}'.format(gene, intronLen))
                    last_exon_end = int(contents[4])
                # mark this as the first gene's exon encountered
                else:
                    last_exon_end = int(contents[4])
    # output the last gene's exon count to stderr
    print('{0}\t{1}'.format(gene, exon_count), file=sys.stderr)
