#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict

def read_bed:
    with open(beda, 'r') as bedfile:
        for bed in bedfile:
            try:
                (chr, start, end, name, score, strand) = bed.strip().split('\t')
            except ValueError:
                'Error reading %s' %(bed)


def main():
    parser = argparse.ArgumentParser(
        description="Given a read alignment file and a base position on reads, "
            "produce a GFF file of where each base pos is located on the reference.")
    parser.add_argument('-b', '--base_pos_info_file', type=str)
    parser.add_argument('-o', '--output_file', type=str)
    parser.add_argument('input_sam_file', type=str)
    args = parser.parse_args()
    queries = load_pos_in_bed_file(args.base_pos_info_file)
    parse_sam(args.input_sam_file, args.output_file, queries)


if __name__ == "__main__":
    main()

