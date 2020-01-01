#!/usr/bin/env python3


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
    addAttr()           Adds new GFF3 attribute
    delAttr()           Delete a GFF3 attribute
    refreshAttrStr()    This needs to be called if changes were made to
                        any of the attributes. It refreshes
    """
    def __init__(self, line=None, **kwargs):
        """GFF3_line is initialized to contain the fields of the GFF3
        line provided as attributes. The attributes are kept in a
        dictionary and a the keys are ordered in a list to preserve
        the order of attributes upon getting the sring representation
        of the line from the GFF3_line object.
        """
        if line == None:
            (self.seqid, 
             self.source, 
             self.type, 
             self.start, 
             self.end, 
             self.score, 
             self.strand, 
             self.phase, 
             self.attributes_str) = [None]*9
            self.attributes_order = []
            self.attributes = {}
        else:
            (self.seqid, 
             self.source, 
             self.type, 
             self.start, 
             self.end, 
             self.score, 
             self.strand, 
             self.phase, 
             self.attributes_str) = line.strip().split('\t')
            self.start = int(self.start)
            self.end = int(self.end)
            assert self.start <= self.end
            self.coords = (self.start, self.end)
            self.length = self.end - self.start + 1
            attributes_list = self.attributes_str.split(';')
            self.attributes_order = [attr.split('=')[0] for attr in 
                                                               attributes_list]
            self.attributes = {attr.split('=')[0]:attr.split('=')[1] for 
                                                       attr in attributes_list}
        self.line_number = None
        if 'line_number' in kwargs:    
            self.line_number = kwargs['line_number']
        # rename the name attribute so it conforms to GFF3 specifications, 
        # where Name is a reserved attribute key. The version of LTRDigest 
        # makes attribute key name
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

    def addAttr(self, attr_key, attr_val, replace=False, attr_pos=0):
        """ds attribute, default is at the start of the list.
        Default behavior is to add attr_val to a growing
        comma-sep list if attr_key exists. Set replace=True to
        replace existing attr_val.
        """
        if attr_key in self.attributes:
            if replace:
                delAttr(attr_key)
                self.attributes[attr_key] = attr_val
                self.attributes_order.insert(attr_pos, attr_key)
                self.refreshAttrStr()
            # grow list
            else:
                self.attributes[attr_key] = '{0},{1}'.format(
                                           self.attributes[attr_key], attr_val)
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
