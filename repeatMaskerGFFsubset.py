#!/usr/bin/env python3

import sys
import re


def help():
    print('''
    Takes RepeatMasker-generated GFF3 on stdin and writes lines from the
    GFF3 file for each of the following categories: copia, gypsy,
    otherltr, nonltr, mudr, hat, and otherdna.
''', file = sys.stderr)
    exit(0)


if '-h' in sys.argv or '-help' in sys.argv:
    help()
for line in sys.stdin:
    if line.startswith('#'):
        continue
    lowercase_name = line.lower()
    if ('copia' in lowercase_name 
      or 'shacop' in lowercase_name):
        with open('copia.gff', 'a') as out:
            out.write(line)
    elif ('gypsy' in lowercase_name 
      or 'ogre' in lowercase_name):
        with open('gypsy.gff', 'a') as out:
            out.write(line)
    elif ('dirs' in lowercase_name
      or 'erv' in lowercase_name
      or 'bel' in lowercase_name
      or 'ltr' in lowercase_name 
      or 'tto1' in lowercase_name 
      or 'tnt1' in lowercase_name 
      or 'tlc1' in lowercase_name 
      or 'tcn1' in lowercase_name 
      or 'gag' in line 
      or 'zf-CCHC' in line 
      or 'DUF4219' in line 
      or 'RVT' in line 
      or 'Asp_protease' in line 
      or 'RVP' in line):
        with open('otherltr.gff', 'a') as out:
            out.write(line)
    elif 'mudr' in lowercase_name:
        with open('mudr.gff', 'a') as out:
            out.write(line)
    elif 'hat' in lowercase_name:
        with open('hat.gff', 'a') as out:
            out.write(line)
    elif ('helitron' in lowercase_name 
      or 'dna' in lowercase_name 
      or re.search('dna\d', lowercase_name)
      or 'mariner' in lowercase_name 
      or 'tc1' in lowercase_name
      or 'harb' in lowercase_name 
      or 'enspm' in lowercase_name 
      or 'cacta' in lowercase_name
      or 'penelope' in lowercase_name 
      or 'polinton' in lowercase_name 
      or 'maverick' in lowercase_name 
      or 'piggybac' in lowercase_name 
      or 'dada' in lowercase_name):
        with open('otherdna.gff', 'a') as out:
            out.write(line)
    elif ('sine' in lowercase_name
      or 'jock' in lowercase_name 
      or 'cr1' in lowercase_name 
      or 'crack' in lowercase_name 
      or 'daphne' in lowercase_name 
      or 'line' in lowercase_name
      or 'l1' in lowercase_name 
      or 'tx1' in lowercase_name 
      or 'rep' in lowercase_name
      or 'rtex' in lowercase_name
      or 'rte' in lowercase_name):
        with open('nonltr.gff', 'a') as out:
            out.write(line)
