#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict
from collections import namedtuple
from itertools import chain


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
    genes = defaultdict(list)
    with open(filename, 'r') as infile:
        for line in infile:
            dat = line.strip().split('\t')
            ginfo = "\t".join(dat[0:6])
            if len(dat) == 13:
                genes[ginfo].append([int(dat[7]),int(dat[8])])
    return genes

def report_fused(genes):
    for gene in genes:
        locs = genes[gene]
        ginfo = gene.split('\t')
        ranges = merge_ranges(locs)
        if(len(ranges)>1):
            for range in ranges:
                print "%s\t%s\t%s\t%d" %(gene, ginfo[0], range[0], range[1])



def main():
    parser = argparse.ArgumentParser(
        description='Description here ')
    parser.add_argument('ibed_filename', type=str)
    args = parser.parse_args()
    genes = read_ibed_file(args.ibed_filename)
    report_fused(genes)

if __name__ == '__main__':
    main()