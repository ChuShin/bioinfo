#!/usr/bin/env python

import argparse
import pysam
import sys
from collections import defaultdict

"""tenx_sscaffold.py:
Given an intersected BED file produced by bedtools-intersectbed,
extract border barcodes with zero coverage check if possible for superscaffolding."""


def is_valid_sequence(start, end, window):
    if (end - start) > 2 * window:
        return True
    else:
        return False


def get_appended_barcodes_in_region(samfile, contig, start, end, min_score):
    barcode_freq = defaultdict(int);
    barcodes = [];
    for alignment in samfile.fetch(contig, start, end):
        aln_tags = dict(alignment.tags)
        [read, barcode] = alignment.query_name.split("_")
        if aln_tags['AS'] > min_score:
            barcode_freq[barcode] += 1;
    for barcode in barcode_freq:
        if barcode_freq[barcode] > 1:
            barcodes.append(barcode)
    return barcodes


def get_barcode_in_contigs(bam_filename, region_filename, window, min_score):
    barcodes = defaultdict(lambda : defaultdict(list))
    contigs = []
    samfile = pysam.AlignmentFile(bam_filename, 'rb')
    with open(region_filename, 'r') as bed:
        for region in bed:
            try:
                (chr, start, end, contig) = region.strip().split('\t')
                contigs.append([contig, int(end) - int(start) + 1])
            except ValueError:
                print 'invalid line: %s' %(region)
                continue
            if is_valid_sequence(int(start), int(end), window):
                p5_start = int(start)
                p5_end = int(start) + window - 1
                p3_start = int(end) - window
                p3_end = int(end)
                barcodes[contig]['p5'] = \
                    get_appended_barcodes_in_region(samfile, contig, p5_start, p5_end, min_score)
                barcodes[contig]['p3'] = \
                    get_appended_barcodes_in_region(samfile, contig, p3_start, p3_end, min_score)
            else:
                print 'skipped: %s is smaller than 2 * %d min_length' %(contig, window)
    samfile.close()
    return barcodes, contigs

def cmp_barcodes(list1, list2):
    return len(set(list1).intersection(list2))

def check_barcode_pairs(barcodes, contigs):
    for idx, contig_info in enumerate(contigs):
        try:
            contig, contig_len = contig_info
            for idx2 in range(idx+1, len(contigs)):
                ncontig, ncontig_len = contigs[idx+1]
                num_links = cmp_barcodes(barcodes[contig]['p3'],barcodes[ncontig]['p5'])
                num_revlinks = cmp_barcodes(barcodes[contig]['p3'],barcodes[ncontig]['p3'])
               if num_links > 5:
                    print '%s join : %s %s %d %d %d' %(contig, ncontig, num_links, num_revlinks, contig_len)
                if num_revlinks > 5:
                    print '%s revjoin : %s %s %d %d %d' %(contig, ncontig, num_links, num_revlinks, contig_len)
                if (num_links <= 5 and num_revlinks <= 5):
                    print '%s %s %d' %(chr, contig, contig_len)
            except IndexError:
                print '%s %s %d' %(contig, contig_len)


def main():
    parser = argparse.ArgumentParser(
        description='Parser for barcoded BAM')
    parser.add_argument('-ms', '--min_score', type=int, default=100)
    parser.add_argument('-w', '--end_window', type=int, default=10000)
    parser.add_argument('bam_filename', type=str)
    parser.add_argument('region_filename', type=str)
    args = parser.parse_args()
    barcodes, contigs = \
        get_barcode_in_contigs(args.bam_filename, args.region_filename, args.end_window, args.min_score)
    check_barcode_pairs(barcodes, contigs)


if __name__ == '__main__':
    main()
