#!/usr/bin/env python

"""
This program take a bed-like file format with the following columns:
ref_chr,ref_start,ref_end,scaffold,sc_start,sc_end
and output the order of scaffolds in ref_chr based on median ref_starts
and relative orientation of scaffold in ref_chr
"""

import argparse
import sys
from numpy import median
from collections import defaultdict

def get_scaffolds_in_chrs(input_file):
    scaffolds =  defaultdict(lambda : defaultdict(list))
    with open(input_file,'r') as bed:
        for line in bed:
            [chr, chr_start, chr_end, sc, sc_start, sc_end] = \
                line.strip().split('\t')
            scaffolds[sc][chr].append([int(chr_start),int(sc_start)])
    return scaffolds


def filter_promicuous_scaffolds(scaffolds):
    chromosomes = defaultdict(lambda: defaultdict(list))
#    results = defaultdict(lambda: defaultdict(int))
    for sc in scaffolds:
        for chr,coords in scaffolds[sc].items():
            print sc,chr,coords
            sc_chr_coords = []
            sc_coords = []
            sc_orient = '+'
            for coord in coords:
                sc_chr_coords.append(coord[0])
                sc_coords.append(coord[1])
            if len(coords) > 1 and sc_coords[0] > sc_coords[-1]:
                sc_orient = '-'
            print sc,chr,len(coords), sc_orient, median(sc_chr_coords)


def main():
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('bed_filename', type=str)
    args = parser.parse_args()
    scaffolds = get_scaffolds_in_chrs(args.bed_filename)
    chromosomes = filter_promicuous_scaffolds(scaffolds)

if __name__ == '__main__':
    main()
