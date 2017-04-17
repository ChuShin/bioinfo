#!/usr/bin/env python

import argparse
import pysam
import sys
from collections import defaultdict

"""tenx_sscaffold.py:
Given an intersected BED file produced by bedtools-intersectbed,
extract border barcodes with zero coverage check if possible for superscaffolding."""


def is_valid_scaffold(start, end, window):
    if (end - start) > 2 * window:
        return True
    else:
        return False


def get_barcodes_in_region(samfile, chr, start, end, min_score):
    barcode_freq = defaultdict(int);
    barcodes = [];
    for alignment in samfile.fetch(chr, start, end):
        aln_tags = dict(alignment.tags)
        if aln_tags['AS'] > min_score and 'MI' in aln_tags:
            barcode_freq[aln_tags['BX']] += 1;
    for barcode in barcode_freq:
        if barcode_freq[barcode] > 1:
            barcodes.append(barcode)
    return barcodes


def get_barcodes_in_scaffold(bam_filename, region_filename, window, min_score):
    barcodes = defaultdict(lambda : (lambda : (lambda : list(int))))
    samfile = pysam.AlignmentFile(bam_filename, 'rb')
    with open(region_filename, 'r') as bed:
        for region in bed:
            try:
                (chr, start, end, scaffold) = region.strip().split('\t')
            except ValueError:
                print 'invalid line: %s' %(region)
                continue
            if is_valid_scaffold(int(start), int(end), window):
                p5_start = int(start)
                p5_end = int(start) + window - 1
                p3_start = int(end) - window
                p3_end = int(end)
                p5_barcodes = get_barcodes_in_region(samfile, chr, p5_start, p5_end, min_score)
                print p5_barcodes
                p3_barcodes = get_barcodes_in_region(samfile, chr, p3_start, p3_end, min_score)
                print p3_barcodes
            else:
                print 'skipped: %s is smaller than 2 * %d min_length' %(scaffold, window)


    samfile.close()



def main():
    parser = argparse.ArgumentParser(
        description='Parser for barcoded BAM')
    parser.add_argument('-mi', '--molecule_only', type=bool, default=1)
    parser.add_argument('-ms', '--min_score', type=int, default=100)
    parser.add_argument('-w', '--end_window', type=int, default=10000)
    parser.add_argument('bam_filename', type=str)
    parser.add_argument('region_filename', type=str)
    args = parser.parse_args()
    get_barcodes_in_scaffold(args.bam_filename, args.region_filename, args.end_window, args.min_score)


if __name__ == '__main__':
    main()
