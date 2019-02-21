#!/usr/bin/env python

import argparse
import pysam
import sys
from collections import defaultdict

"""tenx_sscaffold.py:
Given a BED file, extract border barcodes."""


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
        if alignment.mapping_quality > 20 and \
                aln_tags['AS'] > min_score and 'MI' in aln_tags:
            barcode_freq[aln_tags['BX']] += 1;
    for barcode in barcode_freq:
        if barcode_freq[barcode] > 1:
            barcodes.append(barcode)
    return list(barcodes)


def get_barcodes_in_scaffold(bam_filename, region_filename, window, min_score):
    barcodes = defaultdict(lambda : defaultdict(lambda : defaultdict(list)))
    chromosomes = defaultdict(list)
    samfile = pysam.AlignmentFile(bam_filename, 'rb')
    with open(region_filename, 'r') as bed:
        for region in bed:
            try:
                (chr, start, end, scaffold) = region.strip().split('\t')
                chromosomes[chr].append([scaffold, int(end) - int(start) + 1])
            except ValueError:
                print 'invalid line: %s' %(region)
                continue
            if is_valid_scaffold(int(start), int(end), window):
                p5_start = int(start)
                p5_end = int(start) + window - 1
                p3_start = int(end) - window
                p3_end = int(end)
                p5 = get_barcodes_in_region(samfile, chr, p5_start, p5_end, min_score)
                for barcode in p5:
                    print "%s\tp5\t%d\t%d\t%s\t%s" %(chr, p5_start, p5_end, scaffold, barcode)
                p3 = get_barcodes_in_region(samfile, chr, p3_start, p3_end, min_score)
                for barcode in p3:
                    print "%s\tp3\t%d\t%d\t%s\t%s" %(chr, p3_start, p3_end, scaffold, barcode)
            else:
                print 'skipped: %s is smaller than 2 * %d min_length' %(scaffold, window)
    samfile.close()
    return barcodes, chromosomes

def cmp_barcodes(list1, list2):
    return len(set(list1).intersection(list2))

def check_barcode_pairs(barcodes, chromosomes):
    for chr in sorted(chromosomes):
        for idx, scaffold_info in enumerate(chromosomes[chr]):
            try:
                scaffold, scaffold_len = scaffold_info
                nscaffold, nscaffold_len = chromosomes[chr][idx+1]
                num_links = cmp_barcodes(barcodes[chr][scaffold]['p3'],barcodes[chr][nscaffold]['p5'])
                num_revlinks = cmp_barcodes(barcodes[chr][scaffold]['p3'],barcodes[chr][nscaffold]['p3'])
                if num_links > 5:
                    print '%s join : %s %s %d %d %d' %(chr, scaffold, nscaffold, num_links, num_revlinks, scaffold_len)
                if num_revlinks > 5:
                    print '%s revjoin : %s %s %d %d %d' %(chr, scaffold, nscaffold, num_links, num_revlinks, scaffold_len)
                if (num_links <= 5 and num_revlinks <= 5):
                        print '%s %s %d' %(chr, scaffold, scaffold_len)
            except IndexError:
                print '%s %s %d' %(chr, scaffold, scaffold_len)

def main():
    parser = argparse.ArgumentParser(
        description='Parser for barcoded BAM')
    parser.add_argument('-mi', '--molecule_only', type=bool, default=1)
    parser.add_argument('-ms', '--min_score', type=int, default=100)
    parser.add_argument('-w', '--end_window', type=int, default=30000)
    parser.add_argument('bam_filename', type=str)
    parser.add_argument('region_filename', type=str)
    args = parser.parse_args()
    barcodes, chromosomes = \
        get_barcodes_in_scaffold(args.bam_filename, args.region_filename, args.end_window, args.min_score)
#    check_barcode_pairs(barcodes, chromosomes)


if __name__ == '__main__':
    main()
