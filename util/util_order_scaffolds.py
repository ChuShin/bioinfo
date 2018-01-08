#!/usr/bin/env python

"""
This program take a bed-like file format with the following columns:
ref_chr,ref_start,ref_end,scaffold,sc_start,sc_end
and output the order of scaffolds in ref_chr based on median ref_starts
and relative orientation of scaffold in ref_chr
"""

import argparse
import sys
import numpy
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

    for sc in scaffolds:
        chr_sort = sorted(scaffolds[sc].items(), key=lambda x: len(x[1]),reverse=True)
        chr = chr_sort[0][0]
        coords = chr_sort[0][1]
        num_best = len(coords)

        if(len(chr_sort) > 1):
            num_runner = len(chr_sort[1][1])
            if float(num_runner)/num_best > 0.6:
                #debug message here
                continue
#        print sc,chr,coords
        sc_chr_coords = []
        sc_coords = []
        sc_orient = '+'
        for coord in coords:
            sc_chr_coords.append(coord[0])
            sc_coords.append(coord[1])
        sc_anchor_pos = median(sc_chr_coords)
        if len(coords) > 1:
            corr = numpy.corrcoef(sc_chr_coords,sc_coords)
    #        print corr[0][1]
            if corr[0][1] < 0:
                sc_orient = '-'
        chromosomes[chr][sc].append([sc_anchor_pos, num_best,sc_orient])
    return chromosomes


def main():
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('bed_filename', type=str)
    args = parser.parse_args()
    scaffolds = get_scaffolds_in_chrs(args.bed_filename)
    chromosomes = filter_promicuous_scaffolds(scaffolds)
    for chr in chromosomes:
        sc_sort = sorted(chromosomes[chr].items(), key=lambda x: x[1])
#        print sc_sort
        for sc_order in sc_sort:
            sc = sc_order[0]
            sc_info = chromosomes[chr][sc]
            print "%s\t%s\t%f\t%s\t%s" %(chr, sc,sc_info[0][0],sc_info[0][1],sc_info[0][2])

if __name__ == '__main__':
    main()
