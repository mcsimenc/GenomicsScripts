#!/usr/bin/env python3

import sys


def help():
	print('''
		Usage:
		------------
		blastFilter.py -b <blast_output> \\
                       -floor <min_identity> \\
                       -ceiling <max_identity>

		Description:
		------------
        Takes as input the output of a BLAST program such as blastn, 
        blastp, etc. run with -outfmt 6 or 7 (tabular format) and
        outputs lines if the percent identity of the high scoring pair
        is between the values provided by -floor and -ceiling.

		Options:
		------------
        -b        <path>    blast table
		-floor    <float>   minimum percent identity required for output
                            (default 0)
		-ceiling  <float>   maximum percent identity allowed for output
                            (default 100)
		''')
	sys.exit(0)


# output help information if command line arguments are missing
args = sys.argv
if '-b' not in args or '-h' in args:
	help()
# read command line arguments
else:
	infl = args[ args.index('-b') + 1 ]
if '-floor' in args:
	floor = float(args[ args.index('-floor') + 1 ])
else:
	floor = float(0)
if '-ceiling' in args:
	ceiling = float(args[ args.index('-ceiling') + 1 ])
else:
	ceiling = float(100)
# for each line in the input file read it then check if the percent
# identity is between floor and ceiling. if it is, output the line
with open(infl) as fl:
	for line in fl:
        if not line.startswith('#'):
            contents = line.strip().split('\t')
            id = float(contents[2])
            if id <= ceiling and id >=floor:
                print(line, end='')
