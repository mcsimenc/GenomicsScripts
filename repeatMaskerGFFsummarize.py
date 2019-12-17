#!/usr/bin/env python3

import sys
import re


def help():
    print('''
    Usage:
    ------------
    repeatMaskerGFFsummarize.py <path>

    Description:
    ------------
    Writes three tables with summarized counts and lengths for
    features in a RepeatMasker-derived GFF3: "full", "short", and
    "supershort". The "full" table has the unaltered GFF3 attributes,
    the "short" table has the same categories as "supershort" except
    all features that would fall into the Other category are listed with
    unaltered GFF3 attributes as names. Review of the short table can 
    thus provide indication if another category should be added to this
    script. The "supershort" table outputs summaries for these
    categories, if present:

        Simple
        Unspecified
        Copia
        Gypsy
        hAT
        Helitron
        Mariner
        MuDR
        Harbinger
        Penelope
        Polinton
        PiggyBac
        tRNA
        rRNA
        DADA
        Other_SINE
        EnSpm/CACTA
        RTEX
        RTE
        Jockey
        Other_LINE
        DIRS
        BEL
        Other_LTR_retrotransposon
        L1
        L2
        ERV
        REP
        SOLA
        Satellite
        Microsatellite
        Other_DNA_transposon

    Options:
    ------------
    -db    <path>   Add extra information will be added to the 
                    output table from Repbase DB in EMBL format
    -bed            BED input format (output from bedtools subtract)
    ''')
    sys.exit(0)


# output help if missing command line arguments or -h is present
if len(sys.argv) < 2 or '-h' in sys.argv:
    help()
# regex patterns to match various formats in gff attributes column
target_re_pattern = re.compile('Target "Motif:(.+)"')
# matches repeatmasker- and repeatrunner-derived features
target_re_pattern_repeatmasker = re.compile('Target=(.+?)\s')
# for matching "simple" repeats in repeatmasker results
name_simple_pattern = re.compile('\([ATCG]+\)n')
# for parsing repeatmasker results from a MAKER-derived gff
target_re_pattern_makerlines = re.compile('Name=species:(.+?)\|genus')
# dictionaries to contain summary results
summary_dct_full = {}
summary_dct_short = {}
summary_dct_supershort = {}
# read command line arguments
input_file_name = sys.argv[1]
# read RepBase EMBL-format database
repbase_db = {}
if '-db' in sys.argv:
    with open(sys.argv[sys.argv.index('-db') + 1]) as repbase_db_fl:
        name = ''
        for line in repbase_db_fl:
            if line.startswith('ID'):
                name = line.strip().lstrip('ID ').split()[0]
                repbase_db[name] = ''
            elif line.startswith('DE'):
                repbase_db[name] += line.strip().lstrip('DE ')
            elif line.startswith('KW'):
                repbase_db[name] += line.strip().lstrip('KW ')
# regex pattern to match an unknown format in MAKER-derived gff
bestblasthit_pat = re.compile('bestblasthit_(.+)\.aln_len')
# read input gff and assign each feature to a category in the full,
# short, and supershort dictionaries
with open(input_file_name) as in_fl:
    for line in in_fl:
        # skip commented lines
        if not line.startswith('#'):
            fields = line.split('\t')
            # input is in BED format
            if '-bed' in sys.argv:
                start = int(fields[1])
                end = int(fields[2])
                length = end - start
                # ignore pfam domain lines
                if fields[7] == 'pfam_domain':
                    continue
                # line is MAKER-derived
                if fields[7] == 'match':
                    if 'bestblast' in fields[9]:
                        name = re.search(bestblasthit_pat, fields[9]).group(1)
                    elif 'genus:Simple' in fields[9]:
                        name = 'Simple'
                    else:
                        try:
                            name = re.search(target_re_pattern, 
                                             fields[9]).group(1)
                        # previous regex search failed to match anything
                        except AttributeError:
                            name = re.search(target_re_pattern_repeatmasker, 
                                             fields[9]).group(1)
                    supershort_name = None
                # repeatmasker lines
                else:
                    name = re.search(target_re_pattern, fields[9]).group(1)
                    supershort_name = None
            # input is in gff3 format
            else:
                start = int(fields[3])
                end = int(fields[4])
                length = end - start + 1
                if 'Simple_repeat' in fields[8]:
                    name = 'Simple'
                elif 'bestblast' in fields[8]:
                    name = re.search(bestblasthit_pat, fields[8]).group(1)
                elif '|genus:Unspecified;' in fields[8]:
                    name = 'Unspecified'
                else:
                    # attempt a series of searches with regex patterns
                    try:
                        name = re.search(target_re_pattern, fields[8]).group(1)
                    except AttributeError:
                        try:
                            name = re.search(target_re_pattern_makerlines, 
                                             fields[8]).group(1)
                        except AttributeError:
                            try:
                                name = re.search(target_re_pattern_repeatmasker, 
                                                 fields[8]).group(1)
                            except AttributeError:
                                name = re.search(bestblasthit_pat, 
                                                 fields[8]).group(1)
                supershort_name = None
            # assign categories for the short and supershort tables
            try:
                # -db was used
                repbase_db[name]
                lowercase_name = name.lower()
                repbase_db_entry_lowercase = repbase_db[name].lower()
                if re.search(name_simple_pattern, name):
                    short_name = 'Simple'
                elif name == 'Simple':
                    short_name = 'Simple'
                elif name == 'Unspecified':
                    short_name = 'Unspecified'
                elif ('copia' in lowercase_name 
                        or 'copia' in repbase_db_entry_lowercase 
                        or 'shacop' in lowercase_name):
                    short_name = 'Copia'
                elif ('gypsy' in lowercase_name 
                        or 'gypsy' in repbase_db_entry_lowercase 
                        or 'ogre' in lowercase_name 
                        or 'ogre' in repbase_db_entry_lowercase):
                    short_name = 'Gypsy'
                elif ('helitron' in lowercase_name 
                        or 'helitron' in repbase_db_entry_lowercase):
                    short_name = 'Helitron'
                elif ('mariner' in lowercase_name 
                        or 'mariner' in repbase_db_entry_lowercase 
                        or 'tc1' in lowercase_name):
                    short_name = 'Mariner'
                elif ('mudr' in lowercase_name 
                        or 'mudr' in repbase_db_entry_lowercase 
                        or 'mutator' in repbase_db_entry_lowercase):
                    short_name = 'MuDR'
                elif ('harb' in lowercase_name 
                        or 'harbinger' in repbase_db_entry_lowercase):
                    short_name = 'Harbinger'
                elif 'trna' in lowercase_name:
                    short_name = 'tRNA'
                elif ('rrna' in lowercase_name 
                        or 'rrna' in repbase_db_entry_lowercase):
                    short_name = 'rRNA'
                elif ('snrna' in lowercase_name 
                        or 'snrna' in repbase_db_entry_lowercase):
                    short_name = 'snRNA'
                elif 'sine' in lowercase_name:
                    short_name = 'Other_SINE'
                elif ('cacta' in repbase_db_entry_lowercase 
                        or 'enspm' in lowercase_name 
                        or 'enspm' in repbase_db_entry_lowercase 
                        or 'cacta' in lowercase_name):
                    short_name = 'EnSpm/CACTA'
                elif 'caulimovir' in repbase_db_entry_lowercase:
                    short_name = 'Caulimovirus'
                elif ('microsatellite' in repbase_db_entry_lowercase 
                        or 'minisat' in lowercase_name):
                    short_name = 'Microsatellite'
                elif 'hat' in lowercase_name or 'hAT' in repbase_db[name]:
                    short_name = 'hAT'
                elif 'novosib' in repbase_db_entry_lowercase:
                    short_name = 'Novosib'
                elif ('penelope' in lowercase_name 
                        or 'penelope' in repbase_db_entry_lowercase):
                    short_name = 'Penelope'
                elif ('polinton' in lowercase_name 
                        or 'polinton' in repbase_db_entry_lowercase):
                    short_name = 'Polinton'
                elif ('maverick' in lowercase_name 
                        or 'maverick' in repbase_db_entry_lowercase):
                    short_name = 'Polinton'
                elif ('piggybac' in lowercase_name 
                        or 'piggybac' in repbase_db_entry_lowercase):
                    short_name = 'PiggyBac'
                elif ('RTEX' in repbase_db[name] or 'rtex' in lowercase_name:
                    short_name = 'RTEX'
                elif 'RTE' in repbase_db[name] or 'rte' in lowercase_name:
                    short_name = 'RTE'
                elif ('jock' in lowercase_name 
                        or 'cr1' in lowercase_name 
                        or 'crack' in lowercase_name 
                        or 'daphne' in lowercase_name 
                        or 'jock' in repbase_db_entry_lowercase 
                        or 'crack' in repbase_db_entry_lowercase 
                        or 'daphne' in repbase_db_entry_lowercase 
                        or 'cr1' in repbase_db_entry_lowercase):
                    short_name = 'Jockey'
                elif 'line' in lowercase_name:
                    short_name = 'Other_LINE'
                elif 'dirs' in lowercase_name:
                    short_name = 'DIRS'
                elif ('dada' in lowercase_name 
                        or 'dada' in repbase_db_entry_lowercase):
                    short_name = 'DADA'
                elif 'bel' in lowercase_name:
                    short_name = 'BEL'
                elif ('L1' in repbase_db[name] 
                        or 'l1' in lowercase_name 
                        or 'tx1' in lowercase_name 
                        or 'tx1' in repbase_db[name]):
                    short_name = 'L1'
                elif 'L2' in repbase_db[name] or 'l2' in lowercase_name:
                    short_name = 'L2'
                elif 'erv' in lowercase_name:
                    short_name = 'ERV'
                elif 'rep' in lowercase_name:
                    short_name = 'REP'
                elif 'sola' in lowercase_name or 'sola' in repbase_db[name]:
                    short_name = 'SOLA'
                elif ('satellite' in repbase_db_entry_lowercase 
                        or 'sat' in lowercase_name):
                    short_name = 'Satellite'
                elif 'non-ltr retrotransposon' in repbase_db_entry_lowercase:
                    short_name = 'Other_non-LTR_retrotransposon'
                elif ('dna transpos' in repbase_db_entry_lowercase 
                        or 'dna' in lowercase_name 
                        or re.search('dna\d', lowercase_name)):
                    short_name = 'Other_DNA_transposon'
                elif ('ltr retrotranspos' in repbase_db_entry_lowercase 
                        or 'ltr' in lowercase_name 
                        or 'tto1' in lowercase_name 
                        or 'tnt1' in lowercase_name 
                        or 'tlc1' in lowercase_name 
                        or 'tcn1' in lowercase_name):
                    short_name = 'Other_LTR_retrotransposon'
                elif 'unspecified' in lowercase_name:
                    short_name = 'Unspecified'
                else:
                    short_name = name
                    supershort_name = 'Other'
                if not supershort_name == 'Other':
                    supershort_name = short_name
            # runs if -db was not specified
            except KeyError:
                lowercase_name = name.lower()
                if re.search(name_simple_pattern, name):
                    short_name = 'Simple'
                elif name == 'Simple':
                    short_name = 'Simple'
                elif name == 'Unspecified':
                    short_name = 'Unspecified'
                elif 'copia' in lowercase_name or 'shacop' in lowercase_name:
                    short_name = 'Copia'
                elif 'gypsy' in lowercase_name or 'ogre' in lowercase_name:
                    short_name = 'Gypsy'
                elif 'hat' in lowercase_name:
                    short_name = 'hAT'
                elif 'helitron' in lowercase_name:
                    short_name = 'Helitron'
                elif 'mariner' in lowercase_name or 'tc1' in lowercase_name:
                    short_name = 'Mariner'
                elif 'mudr' in lowercase_name:
                    short_name = 'MuDR'
                elif 'harb' in lowercase_name:
                    short_name = 'Harbinger'
                elif 'penelope' in lowercase_name:
                    short_name = 'Penelope'
                elif 'polinton' in lowercase_name:
                    short_name = 'Polinton'
                elif 'piggybac' in lowercase_name:
                    short_name = 'PiggyBac'
                elif 'trna' in lowercase_name:
                    short_name = 'tRNA'
                elif 'rrna' in lowercase_name:
                    short_name = 'rRNA'
                elif 'dada' in lowercase_name:
                    short_name = 'DADA'
                elif 'sine' in lowercase_name:
                    short_name = 'Other_SINE'
                elif 'enspm' in lowercase_name or 'cacta' in lowercase_name:
                    short_name = 'EnSpm/CACTA'
                elif 'rtex' in lowercase_name:
                    short_name = 'RTEX'
                elif 'rte' in lowercase_name:
                    short_name = 'RTE'
                elif ('jock' in lowercase_name 
                        or 'cr1' in lowercase_name 
                        or 'crack' in lowercase_name 
                        or 'daphne' in lowercase_name):
                    short_name = 'Jockey'
                elif 'line' in lowercase_name:
                    short_name = 'Other_LINE'
                elif 'dirs' in lowercase_name:
                    short_name = 'DIRS'
                elif 'bel' in lowercase_name:
                    short_name = 'BEL'
                elif ('ltr' in lowercase_name 
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
                    short_name = 'Other_LTR_retrotransposon'
                elif 'l1' in lowercase_name or 'tx1' in lowercase_name:
                    short_name = 'L1'
                elif 'l2' in lowercase_name:
                    short_name = 'L2'
                elif 'erv' in lowercase_name:
                    short_name = 'ERV'
                elif 'rep' in lowercase_name:
                    short_name = 'REP'
                elif 'sola' in lowercase_name:
                    short_name = 'SOLA'
                elif 'sat' in lowercase_name:
                    short_name = 'Satellite'
                elif 'minisat' in lowercase_name:
                    short_name = 'Microsatellite'
                elif ('dna' in lowercase_name 
                        or re.search('dna\d', lowercase_name)):
                    short_name = 'Other_DNA_transposon'
                elif 'unspecified' in lowercase_name:
                    short_name = 'Unspecified'
                else:
                    short_name = name
                    supershort_name = 'Other'
                if not supershort_name == 'Other':
                    supershort_name = short_name
            # for full output
            if name in summary_dct_full:
                summary_dct_full[name]['count'] += 1
                summary_dct_full[name]['length'] += length
            else:
                summary_dct_full[name] = {'count':1, 'length':length}
            # for short output
            if short_name in summary_dct_short:
                summary_dct_short[short_name]['count'] += 1
                summary_dct_short[short_name]['length'] += length
            else:
                summary_dct_short[short_name] = {'count':1, 'length':length}
            # for supershort output
            if supershort_name in summary_dct_supershort:
                summary_dct_supershort[supershort_name]['count'] += 1
                summary_dct_supershort[supershort_name]['length'] += length
            else:
                summary_dct_supershort[supershort_name] = {
                                                    'count':1, 'length':length}
# write full table
with open('{0}.summary_full'.format(input_file_name), 'w') as out_fl_full:
        if '-db' in sys.argv:
            out_fl_full.write('Name\tCount\tLength\tDescription\n')
            for name in summary_dct_full:
                try:
                    out_fl_full.write('{0}\t{1}\t{2}\t{3}\n'.format(
                                               name, 
                                               summary_dct_full[name]['count'], 
                                               summary_dct_full[name]['length'], 
                                               repbase_db[name]))
                except KeyError:
                    out_fl_full.write('{0}\t{1}\t{2}\t{3}\n'.format(
                                               name, 
                                               summary_dct_full[name]['count'], 
                                               summary_dct_full[name]['length'], 
                                               ('No description found. Name '
                                                'may not exactly match ID in '
                                                'rebpase EMBL')))
        else:
            # default, if -db is not specified
            out_fl_full.write('Name\tCount\tLength\n')
            for name in summary_dct_full:
                out_fl_full.write('{0}\t{1}\t{2}\n'.format(name, 
                                            summary_dct_full[name]['count'], 
                                            summary_dct_full[name]['length']))
# output short table
with open('{0}.summary_short'.format(input_file_name), 'w') as out_fl_short:
        if '-db' in sys.argv:
            out_fl_short.write('Name\tCount\tLength\tDescription\n')
            for name in summary_dct_short:
                try:
                    out_fl_short.write('{0}\t{1}\t{2}\t{3}\n'.format(name, 
                                            summary_dct_short[name]['count'], 
                                            summary_dct_short[name]['length'], 
                                            repbase_db[name]))
                except KeyError:
                    out_fl_short.write('{0}\t{1}\t{2}\t{3}\n'.format(name, 
                                            summary_dct_short[name]['count'], 
                                            summary_dct_short[name]['length'], 
                                            ('No description found. Name may '
                                             'not exactly match ID in rebpase '
                                             'EMBL')))
        else:
            # default, if -db is not specified
            out_fl_short.write('Name\tCount\tLength\n')
            for name in summary_dct_short:
                out_fl_short.write('{0}\t{1}\t{2}\n'.format(name, 
                                            summary_dct_short[name]['count'], 
                                            summary_dct_short[name]['length']))
# output supershort table
with open('{0}.summary_supershort'.format(
                                   input_file_name), 'w') as out_fl_supershort:
        if '-db' in sys.argv:
            out_fl_supershort.write('Name\tCount\tLength\tDescription\n')
            for name in summary_dct_supershort:
                try:
                    out_fl_supershort.write('{0}\t{1}\t{2}\t{3}\n'.format(name, 
                                        summary_dct_supershort[name]['count'], 
                                        summary_dct_supershort[name]['length'], 
                                        repbase_db[name]))
                except KeyError:
                    out_fl_supershort.write('{0}\t{1}\t{2}\t{3}\n'.format(name, 
                                        summary_dct_supershort[name]['count'], 
                                        summary_dct_supershort[name]['length'], 
                                        ('No description found. Name may not '
                                        'exactly match ID in rebpase EMBL')))
        else:
            # default, if -db is not specified
            out_fl_supershort.write('Name\tCount\tLength\n')
            for name in summary_dct_supershort:
                out_fl_supershort.write('{0}\t{1}\t{2}\n'.format(name, 
                                       summary_dct_supershort[name]['count'], 
                                       summary_dct_supershort[name]['length']))
