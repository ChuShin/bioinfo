#!/usr/bin/env python


import argparse
import sys
import argparse
import sys
from collections import defaultdict



Ibed = namedtuple('ibed', 'chr, gs, ge, gid, gscore, gstrand, '
                                    'mchr, mgs, mge, mgid, mgscore, mgstrand,'
                                    'ovl')
flatten = chain.from_iterable


def merge_ranges(data):
    data = sorted(flatten(((start, 1), (stop, -1)) for start, stop in data))
    merged = []
    start, end, count = 0, 0, 0
    for pos, side in data:
        if count == 0:
            start = pos
        count += side
        if count == 0:
            end = pos
            merged.append([start, end])
    return merged

def read_ibed_file(filename):
    genes = defaultdict(lambda: defaultdict(int))

    with open(filename, 'r') as infile:
        for line in infile:
            dat = line.strip().split('\t')
            if len(dat) == 13:
                alns.append(Ibed._make(dat))
    return alns




def main():
    parser = argparse.ArgumentParser(
        description='Description here ')
    parser.add_argument('ibed_filename', type=str)
    args = parser.parse_args()
    genes = read_ibed_file(args.ibed_filename)
    lift_over(args.bed_filename,chrs)

if __name__ == '__main__':
    main()