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
        [read, barcode] = alignment.query_name.split("_")
        if aln_tags['AS'] > min_score and 'MI' in aln_tags:
            barcode_freq[aln_tags['BX']] += 1;
    for barcode in barcode_freq:
        if barcode_freq[barcode] > 1:
            barcodes.append(barcode)
    return barcodes


def get_appended_barcodes_in_region(samfile, chr, start, end, min_score):
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
    scaff_barcodes = defaultdict(lambda : defaultdict(lambda : defaultdict(list)))
    scaff_lengths = []
    samfile = pysam.AlignmentFile(bam_filename, 'rb')
    with open(region_filename, 'r') as bed:
        for region in bed:
            try:
                (chr, start, end, scaffold) = region.strip().split('\t')
                scaff_lengths.append([scaffold, int(end) - int(start) + 1])
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
                scaff_barcodes[scaffold]['p5']= p5
                scaff_barcodes[scaffold]['p3'] = p3
            else:
                print 'skipped: %s is smaller than 2 * %d min_length' %(scaffold, window)
    samfile.close()
    return scaff_barcodes, scaff_lengths


def cmp_barcodes(list1, list2):
    return len(set(list1).intersection(list2))

def check_barcode_pairs(barcodes, scaffs):
    for sc1_ind in range(0,len(scaffs)-1):
        for sc2_ind in range(sc1_ind+1,len(scaffs)):
            sc1, sc1_length = scaffs[sc1_ind]
            sc2, sc2_length = scaffs[sc2_ind]
            num_links = cmp_barcodes(barcodes[sc1]['p3'], barcodes[sc2]['p5']) \
                        + cmp_barcodes(barcodes[sc2]['p3'], barcodes[sc1]['p5'])
            num_revlinks = cmp_barcodes(barcodes[sc1]['p5'], barcodes[sc2]['p5']) \
                        + cmp_barcodes(barcodes[sc1]['p3'], barcodes[sc2]['p3'])
            if num_links > 5 or num_revlinks > 5:
                print '%s\t%d\t%s\t%d\t%d\t%d' %(sc1, sc1_length, sc2, sc2_length, num_links, num_revlinks)


def main():
    parser = argparse.ArgumentParser(
        description='Parser for barcoded BAM')
    parser.add_argument('-mi', '--molecule_only', type=bool, default=1)
    parser.add_argument('-ms', '--min_score', type=int, default=100)
    parser.add_argument('-w', '--end_window', type=int, default=30000)
    parser.add_argument('bam_filename', type=str)
    parser.add_argument('region_filename', type=str)
    args = parser.parse_args()

    scaff_barcodes, scaff_lengths = \
        get_barcodes_in_scaffold(args.bam_filename, args.region_filename, args.end_window, args.min_score)
    check_barcode_pairs(scaff_barcodes, scaff_lengths)


if __name__ == '__main__':
    main()
