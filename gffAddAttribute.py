#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    gffAddAttr.py -gff <filepath> -attr <string> -map <filepath> \\
                      -mapKey <string> > output.gff

    Description: 
    ------------
    Adds a new attribute key-value pair to field 9 of the input GFF3
    file. Four arguments are required: -gff, specifying the input GFF3,
    -attr, specifying the name of the key of the new attribute, -mapKey,
    referring to an attribute key present in the original GFF3, and
    -map, a map of values of the key specified by -mapKey, such as
    a list of gene names, to the value of the new attribute keys. The
    purpose of the -map file is to provide a way to map the new
    attribute values to specific features in the input GFF3.

    e.g.
        -mapKey ID -attr function

    original GFF3 line attributes:
    ------------------------------
        ID=Gene1
        ID=Gene2

    the file given by -map:
    ------------------------------
        Gene1   protease
        Gene2   intracellular_transport

    attributes of the modified GFF3:
    ------------------------------
        ID=Gene1;function=protease
        ID=Gene2;function=intracellular_transport
    

    Options:
    ------------
    -attr    <string>    New attribute key name. Must be unique among
                         existing attribute keys in the GFF3 lines.

    -gff    <filepath>    The GFF3 file to modify

    -map    <filepath>    A two-column tab-delimited file where first 
                          column contains the values of the attribute 
                          (in the original GFF3) specified by -mapKey 
                          and the second column to the value for the 
                          new attribute key.

    -mapKey    <string>   An attribute key present in the original GFF3
                          file, the values of which are present in the
                          first column of the file specified by -map

    -restrictType <string> Only add new attribute to features of this 
                           type

    -replace              Replace existing attribute if present

    -replaceIfNone        Replace existing attribute only if present
                          and its value is 'None'

    -v                    Verbose reporting on stderr
        ''', file=sys.stderr)
    sys.exit(0)


class GFF3_line:
    """
    Attributes:
    ------------
    field0, ... , field8 - string
    attributes - dictionary

    Methods:    
    ------------
    str()               Prints GFF3 line.
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

    def refreshAttrStr(self, attrOrder=None):
        """If the attributes dictionary or attributes_order has been 
        altered this should be called to update attributes_str.
        """
        self.attributes_str = ';'.join(['='.join(
             [attr, self.attributes[attr]]) for attr in self.attributes_order])


# print help information if not eanough arguments are present
args =sys.argv
if len(sys.argv) < 9 or '-h' in args:
    help()
attrMap = {}
with open(sys.argv[sys.argv.index('-map') + 1]) as mapFl:
    for line in mapFl:
        try:
            gene, annot = line.strip().split('\t')
            attrMap[gene] = annot
        except:
            continue
mapKey = sys.argv[sys.argv.index('-mapKey')+1]
newAttrKey = sys.argv[sys.argv.index('-attr')+1]
gffFilepath = sys.argv[sys.argv.index('-gff')+1]
if '-restrictType' in args:
    restrictType = sys.argv[sys.argv.index('-restrictType')+1]
with open(gffFilepath) as gffFl:
    for line in gffFl:
        if line.startswith('#'):
            print(line, end='')
        else:
            gffLine = GFF3_line(line)
            if '-restrictType' in args:
                if gffLine.type != restrictType:
                    print(line, end='')
                    continue
            if newAttrKey in gffLine.attributes:
                if (not '-replace' in args 
                    and ('-replaceIfNone' in args 
                    and gffLine.attributes[newAttrKey]) != 'None'):
                    if '-v' in args:
                        print(('{0}\nAbove line in GFF3 input already has '
                              'attribute {1}. Continuing').format(str(gffLine), 
                                                  newAttrKey), file=sys.stderr)
                    print(line, end='')
                    continue
            try:
                featureIdentifier = gffLine.attributes[mapKey]
            except KeyError:
                if '-v' in args:
                    print('{0}\n-mapKey {1} not in above GFF3 line.'.format(
                                        str(gffLine), mapKey), file=sys.stderr)
                    print('--------', file=sys.stderr)
                newAttrValue = "None"
                gffLine.attributes[newAttrKey] = newAttrValue
                if newAttrKey not in gffLine.attributes_order:
                    gffLine.attributes_order.append(newAttrKey)
                try:
                    gffLine.refreshAttrStr()
                except:
                    print(gffLine.attributes, file=sys.stderr)
                    print(gffLine.attributes_order, file=sys.stderr)
                    sys.exit()
                print(str(gffLine))
                continue
            if featureIdentifier in attrMap:
                newAttrValue = attrMap[featureIdentifier]
                gffLine.attributes[newAttrKey] = newAttrValue
                if newAttrKey not in gffLine.attributes_order:
                    gffLine.attributes_order.append(newAttrKey)
                gffLine.refreshAttrStr()
                print(str(gffLine))
            else:
                if '-v' in args:
                    print(('{0}\nNo item in -map matching {1} from the above '
                                            'GFF3 line').format(str(gffLine), 
                                                           featureIdentifier), 
                                                               file=sys.stderr)
                    print('--------', file=sys.stderr)
                newAttrValue = "None"
                gffLine.attributes[newAttrKey] = newAttrValue
                if newAttrKey not in gffLine.attributes_order:
                    gffLine.attributes_order.append(newAttrKey)
                gffLine.refreshAttrStr()
                print(str(gffLine))
