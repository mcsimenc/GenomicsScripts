#!/usr/bin/env python3

import sys
import statistics


def help():
    print('''
    Usage:
    ------------
    meanMedianMinMax < <path>

    Description:
    ------------
    Calculates min, max, mean, and median from a list of numbers 
    given on stdin and outputs to stdout.
    ''')
    sys.exit(0)


if '-h' in sys.argv:
    help()
values = []
for line in sys.stdin:
    values.append(float(line.strip()))
print('min\t{0}'.format(min(values)))
print('max\t{0}'.format(max(values)))
print('mean\t{0}'.format(statistics.mean(values)))
print('median\t{0}'.format(statistics.median(values)))
print('total\t{0}'.format(sum(values)))
