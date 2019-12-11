#!/usr/bin/env python3

import sys

def help():
    print('''
    Usage:
    ------------
    gffv2Exonerate2gff3.py < input.gff2 > output.gff3
    
    Description:
    ------------
    Exonerate's est2genome, protein2genome outputs GFF as GFF2. This
    script takes an Exonerate GFF2 file as input and outputs a GFF3
    format file.
    ''')
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

	def addAttr(self, attr_key, attr_val, replace=False, attr_pos=0):
		"""
		Add key-value pair to the GFF3 attribute column, at the
        beginning of the list by default. If the key exists then the
        value is concatenated with the current value separated by a
        comma, unless replace=True in which case the existing value is
        replaced.
		"""
		if attr_key in self.attributes:
            # replace the current value of this key with then new value
			if replace:
				delAttr(attr_key)
				self.attributes[attr_key] = attr_val
				self.attributes_order.insert(attr_pos, attr_key)
				self.refreshAttrStr()
            # this key already exists in the attributes of the GFF3 file
            # add this value to the existing key's values with commas
            # separating values
			else:
				self.attributes[attr_key] = '{0},{1}'.format(
                                                     self.attributes[attr_key], 
                                                     attr_val)
				self.refreshAttrStr()
		else:
			self.attributes[attr_key] = attr_val
			self.attributes_order.insert(attr_pos, attr_key)
			self.refreshAttrStr()

	def delAttr(self, attr_key):
		"""Deletes attribute."""
		del self.attributes[attr_key]
		self.attributes_order.pop(self.attributes_order.index(attr_key))
		self.refreshAttrStr()

    def refreshAttrStr(self):
        """If the attributes dictionary or attributes_order has been 
        altered this should be called to update attributes_str.
        """
        self.attributes_str = ';'.join(['='.join(
             [attr, self.attributes[attr]]) for attr in self.attributes_order])


def convert_GFF2_to_GFF3(line):
    """
    Takes a GFF2-format line as input and returns a GFF3_line object
    """
    gff3 = GFF3_line()
    if len(line.strip().split('\t')) == 9:
        gff3.seqid, 
        gff3.source, 
        gff3.type, 
        gff3.start, 
        gff3.end, 
        gff3.score, 
        gff3.strand, 
        gff3.phase, 
        attr = line.strip().split('\t')
        if gff3.type == 'similarity':
            return None
        attr = attr.split(';')
        for pair in attr:
            k,v = pair.split()
            gff3.attributes[k] = v
            gff3.attributes_order.append(k)
            gff3.refreshAttrStr()
    elif len(line.strip().split('\t')) == 8:
        gff3.seqid, 
        gff3.source, 
        gff3.type, 
        gff3.start, 
        gff3.end, 
        gff3.score, 
        gff3.strand, 
        gff3.phase = line.strip().split('\t')
        gff3.attributes['ID'] = '.'
        gff3.attributes_order.append('ID')
        if gff3.type == 'similarity':
            return None

    return gff3
    

if __name__ == '__main__':
    args = sys.argv
    # output help information if asked for
    if '-h' in args:
        help()
    begin = False
    print('##gff-version 3')
    feats = {}
    for line in sys.stdin:
        # before the gff2 feature lines begin there are non-feature
        # lines that are not commented. The first commented lines are
        # immediately preceding the feature lines. Ignore other
        # commented lines
        if line.startswith('#'):
            begin = True
            continue
        elif begin:
            gff3 = convert_GFF2_to_GFF3(line)
            if not gff3 == None:
                scaf = gff3.seqid
                coords = (gff3.start,gff3.end)
                # Combine feats with identical type and coords, adding
                # evidence (transcript) to master record
                if scaf in feats:
                    if coords in feats[scaf]:
                        if gff3.type in feats[scaf][coords]:
                            if 'sequence' in gff3.attributes:
                                feats[scaf][coords][gff3.type].addAttr(
                                          attr_key='sequence', 
                                          attr_val=gff3.attributes['sequence'])
                        else:
                            feats[scaf][coords][gff3.type] = gff3
                    else:
                        feats[scaf][coords] = { gff3.type:gff3 }
                else:
                    feats[scaf] = {coords:{gff3.type:gff3}}
    # print GFF3 lines
    for scaf in feats:
        for t in feats[scaf]:
            for c in feats[scaf][t]:
                print(str(feats[scaf][t][c]))
