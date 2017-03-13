#!/usr/bin/env python

import argparse
import sys
from collections import namedtuple
from collections import defaultdict

"""calculate_percent_methylation.py:
Given an intersected BED file produced by bedtools-intersectbed,
calculate the methylation percentage for each feature (column4)."""


Methyl_Data = namedtuple('Methyl_Data', 'chr, start, end, feature, score, '
                                    'strand, bs_chr, bs_start, bs_end, '
                                    'bs_context, bs_type, bs_strand, overlap')


def load(intersect_file):
    pass


def main():
    parser = argparse.ArgumentParser(
        description="Given an intersect bed output file, produce a summary of "
                    "%methyl for each input feature.")
    parser.add_argument('-m', '--min_obs_nuc', type=int, default=10)
    parser.add_argument('-o', '--output_file', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('intersect_file', type=str)
    args = parser.parse_args()
    intersect_features = load(args.intersect_file)
    summarize_ref_simple(alns, refs,
                     args.min_obs_nuc)

if __name__ == "__main__":
    main()