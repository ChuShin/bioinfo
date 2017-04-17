#!/usr/bin/env python

import argparse
import pysam
import sys
from collections import defaultdict

"""tenx_sscaffold.py:
Given an intersected BED file produced by bedtools-intersectbed,
extract border barcodes with zero coverage check if possible for superscaffolding."""


def read_bam(bam_filename, region_filename, min_score):
    samfile = pysam.AlignmentFile(filename, 'rb')
    with open(region_filename, 'r') as bed:
        for region in bed:
            try:
                (chr, start, end) = region.strip().split('\t')
            except ValueError:
                continue
            alignments = samfile.fetch(chr, start, end)
            print alignments


    samfile.close()



def main():
    parser = argparse.ArgumentParser(
        description='Parser for barcoded BAM')
    parser.add_argument('-mi', '--molecule_only', type=bool, default=1)
    parser.add_argument('-ms', '--min_score', type=int, default=100)
    parser.add_argument('bam_filename', type=str)
    parser.add_argument('region_filename', type=str)
    args = parser.parse_args()
    read_bam(args.bam_filename, args.region_filename, args.min_score)


if __name__ == '__main__':
    main()