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
        description="")
    parser.add_argument('-a', '--abed_file', type=str)
    parser.add_argument('-b', '--bbed_file', type=str)
    args = parser.parse_args()
    read_bed(args.abed_file)


if __name__ == "__main__":
    main()

