#!/usr/bin/env python3
import re
import os
import sys
import subprocess
from Bio import SeqIO


def help():
    print('''
    Usage:
    ------------
    gff2fasta.py -gff <path> -fasta <path> [-attr <str>]
    
    Description:
    ------------
    BEDtools must be in the PATH and Biopython must be installed.
    Extracts sequences from the FASTA file provided by -fasta for 
    features in the GFF3 file provided by -gff. The ninth column in 
    GFF3 files contains attributes such as ID and Name. By default, the
    value for the ID attribute is used as the sequence name in the 
    output FASTA file. An alternative attribute key from which to 
    derive the name can be specified using -attr. The output FASTA is 
    named the same as the input GFF3 with the added suffix .fasta. By 
    default the sequences in the FASTA are treated as the + strand; 
    - strand features from the GFF3 are reverse complemented.
    
    Options:
    ------------
    -attr <str>    Attribute from which to derive FASTA sequence name
        ''')
    sys.exit(0)


def bedtoolsid2attr(gff_flpath, 
                    attr='ID', 
                    strand=False, 
                    lstrip=None):
    """Creates a map of bedtools getfasta-stype sequence names to new
    sequence names obtained from the specified attribute from the GFF3
    file and returns it as a dictionary.
    """
    attr_pat = re.compile('{0}=(.+?)(;|$)'.format(attr))
    map_dct = {}
    with open(gff_flpath) as in_fl:
        for line in in_fl:
            # ignore commented lines
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                seqid = fields[0]
                start = str(int(fields[3]) - 1)
                end = fields[4]
                strand = fields[6]
                attr_value = re.search(attr_pat, fields[8]).group(1)
                # remove characters from the left side of the string
                if not lstrip == None:
                    attr_value = attr_value.lstrip(lstrip)
                if strand:
                    # the bedtools getfasta output sequence names are 
                    # different if -s is used (strand-aware extraction)
                    map_dct['{0}:{1}-{2}({3})'.format(seqid, start, end, 
                                                          strand)] = attr_value
                # strand-unaware extraction was performed
                else:
                    map_dct['{0}:{1}-{2}'.format(seqid, start, 
                                                             end)] = attr_value
    # return a dictionary containing the map of old sequence ids to the 
    # to new sequence ids
    return map_dct

def rename_fasta_seq_headers(in_flpath, in_header_pattern, header_to_name_map, out_flpath):
    """
    I/O for changing fasta sequence headers. Uses regex matching and dictionary associating current with new name.
    in_fasta_filepath - the fasta file with headers to be renamed.
    in_header_pattern - a regex that will match the part of the header that corresponds to keys in header_to_name_map
    header_to_name_map - a dict with keys as part of in_fasta_filepath seq headers that match in_header_pattern and values as the abbreviated LTR-RT ID (e.g. 24_1, which corresponds to LTR_retrotransposon24)

    This function depends on the re and Bio.SeqIO modules
    """
    fasta_file = list(SeqIO.parse(in_flpath, 'fasta'))
    in_header_pattern = re.compile(in_header_pattern)

    for seq in fasta_file:
        seq_id = re.search(in_header_pattern, seq.id).group(1)
        seq.id = header_to_name_map[seq_id]
        seq.description = ''
    
    SeqIO.write(fasta_file, out_flpath, 'fasta')

def makecall(call, stdout=None, stderr=None, stdin=None):
    """Simplifies running system calls using the subprocess module. 
    stdout and sterr are written to files automatically.
    """
    # streams are ignored
    if stdout == None and stderr == None and stdin == None:
        subprocess.call(call)
    elif stdout != None:
        with open(stdout, 'w') as outfl:
            if stderr != None:
                with open(stderr, 'w') as errfl:
                    # stderr is written to a new file
                    if stdin == None:
                        subprocess.call(call, stdout=outfl, stderr=errfl)
                    # receives a stream from stdin and writes stdin and 
                    # stdout are written to files
                    else:
                        with open(stdin, 'r') as inFl:
                            subprocess.call(call, stdin=inFl, stdout=outfl, 
                                                                  stderr=errfl)
            elif stderr == None:
                # stdout is written to a file
                if stdin == None:
                    subprocess.call(call, stdout=outfl)
                # receives a stream from stdin and writes stdout to a 
                # new file
                else:
                    with open(stdin, 'r') as inFl:
                        subprocess.call(call, stdin=inFl, stdout=outfl)

    elif stderr != None and stdout == None:
        with open(stderr, 'w') as errfl:
            # stderr is written to a file
            if stdin == None:
                subprocess.call(call, stderr=errfl)
            else:
                # receives a stream from stdin and writes stdout to a 
                # new file
                with open(stdin, 'r') as inFl:
                    subprocess.call(call, stdin=inFl,  stderr=errfl)
    # receives a stream  from stdin
    elif stdin != None and stderr == None and stdout == None:
        with open(stdin, 'r') as inFl:
            subprocess.call(call, stdin=inFl)

def ChangeFastaHeaders(inputFastaPath, inputGFFpath, attribute='ID'):
    """Creates map of bedtools getfasta-style sequence names to the new
    names obtained from the input GFF3. Then a new FASTA file is 
    written.
    """
    # create a map of sequence ids
    bedtoolsIDmap = bedtoolsid2attr(inputGFFpath, attr=attribute)
    # write fasta with new sequence headers obtained fro the map
    newFasta = '{0}.new.tmp'.format(inputFastaPath)
    header_pattern='(.+?:\d+?-\d+?)(?:$|\D)'
    rename_fasta_seq_headers(inputFastaPath, header_pattern, bedtoolsIDmap, 
                                                                      newFasta)
    os.replace(newFasta, inputFastaPath)


def removeRedundant(fastaPth):
    """Removes all but the first sequence with the same name as other
    sequences in the FASTA file. This can happen if multiple lines in
    the GFF3 input have the same sequence but it can also happen if
    more than one sequence in the GFF3 have the same value for the
    attribute used to obtain sequence names from (ID by default).
    If sequences are removed here they are reported.
    """
    seqnames = set()
    tmp = '{0}.nonredundant'.format(fastaPth)
    if os.path.isfile(tmp):
        os.remove(tmp)
    with open(tmp, 'w') as outFl:
        with open(fastaPth, 'r') as inFl:
            skip = False
            for line in inFl:
                # skip commented lines
                if line.startswith('>'):
                    if line in seqnames:
                        print(('Sequence not output due to '
                                          'redundant name:\t{}').format(line),
                               file=sys.stderr)
                        skip = True
                        continue
                    else:
                        skip = False
                        seqnames.add(line)
                        outFl.write(line)
                else:
                    if skip:
                        continue
                    else:
                        outFl.write(line)
    os.rename(tmp, fastaPth)


def getfasta(inGFFpth, fastaRefPth, outFastaPth, headerKey='ID'):
    """Runs bedtools getfasta to extract sequences from the input
    FASTA file. Default bedtools names are renamed according to the
    attribute specified.
    """
    call = ['bedtools', 'getfasta', 
            '-fi', fastaRefPth, 
            '-s', 
            '-bed', inGFFpth ]
    # run bedtools getfasta
    try:
        makecall(call, stdout=outFastaPth)
    except FileNotFoundError:
        print(("BEDtools was not found. Make sure you have bedtools in your "
                                            "PATH."), file=sys.stderr)
        exit(1)
    # rename fasta sequences
    ChangeFastaHeaders(outFastaPth, inGFFpth, attribute=headerKey)
    # remove redundantly named sequences from the fasta file
    removeRedundant(outFastaPth)

if __name__ == '__main__':
    # Parse args
    args = sys.argv
    if "-help" in args or "-h" in args or len(args) < 5:
        help()
    if '-attr' in args:
        attr = args[args.index('-attr')+1]
    else:
        attr = 'ID'
    fastaPth = args[args.index('-fasta')+1]
    gffPth = args[args.index('-gff')+1]
    # Get sequences
    getfasta(gffPth, fastaPth, '{0}.fasta'.format(gffPth), headerKey=attr)
