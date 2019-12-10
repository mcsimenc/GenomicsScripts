#!/usr/bin/env python3

import sys


def help():
    print('''
     Usage:
     ---------
     coverage2circostrack.py  <genomecov_output>  > circos.linetrack
     
     Input:
     ---------
     The output of bedtools genomecov -ibam <bam> -d:
    
     chr    pos    doc
    
     Output:
     ---------
     A line track for Circos. If the final bin is smaller than the
     specified bin size then average depth of coverage may appear more
     extreme than is estimated for full-size bins.
    
     chr     start    end    average_depth_of_coverage
    
     Options:
     ---------
     -bin   <int>   bp in each bin. Default (100000)
 ''')


if __name__ == '__main__':
    # output help information if not enough command line arguments are
    # specified
    args = sys.argv
    if len(args) < 2:
        help()
        sys.exit()
    # read command line arguments
    if '-bin' in args:
        binsize = int(args[args.index('-bin') + 1])
    else:
        binsize = 100000
    with open(args[1],'r') as inFl:
        ct = 0
        scaf_pos = 1
        sum_depth = 0
        last_scaf = None
        for line in inFl:
            ct += 1
            scaf, pos, depth = line.strip().split()
            if last_scaf == None:
                last_scaf = scaf
            # end of a bin. calculate average depth of coverage and
            # output line for circos line track
            if ct == binsize:
                avg_depth = sum_depth / binsize
                print('{0}\t{1}\t{2}\t{3:.4f}'.format(last_scaf, 
                                                      scaf_pos - binsize, 
                                                      scaf_pos - 1, 
                                                      avg_depth))
                ct = 0
                avg_depth = 0
                sum_depth = 0
            # end of a scaffold. calculate average depth of coverage for
            # the last possibly shorter bin and output circos line track
            # line
            if scaf != last_scaf:
                avg_depth = sum_depth / ((scaf_pos) - (scaf_pos - ct - 1))
                print('{0}\t{1}\t{2}\t{3:.4f}'.format(last_scaf, 
                                                      scaf_pos-ct, 
                                                      scaf_pos, 
                                                      avg_depth))
                last_scaf = scaf
                scaf_pos = 0
                avg_depth = 0
                ct = 0
                sum_depth = 0
            scaf_pos += 1
            depth = int(depth)
            sum_depth += depth
        # last bin of the last scaffold in the file
        avg_depth = sum_depth / ((scaf_pos - 1) - (scaf_pos - 1 - ct - 1))
        print('{0}\t{1}\t{2}\t{3:.4f}'.format(last_scaf, 
                                              scaf_pos - 1 - ct - 1, 
                                              scaf_pos - 1, 
                                              avg_depth))
