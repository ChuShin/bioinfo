#!/usr/bin/env python

import argparse
import sys
import gzip
from collections import defaultdict

"""util_check_index.py:
Given a FASTQ file print the total number of read count and top N most
abundant indices."""


def print_indices_summary(fastq_file, num_indices):
    read_count = 0
    rank = 1
    index_lib = defaultdict(int)
    with gzip.open(fastq_file, 'r') as infile:
        for line in infile:
            if line.startswith('@'):
                [header, index] = line.strip().split(' ')
                read_count += 1
                index_lib[index] += 1
    for index in sorted(index_lib, key=index_lib.get, reverse=True)[0:num_indices]:
        print '%s\t%\t%s\t%d' %(fastq_file, rank, index, index_lib[index])
        rank += 1
    print '%s\t%d' %(fastq_file, read_count)


def main():
    parser = argparse.ArgumentParser(
        description='Given a FASTQ file print total read count and top N '
                    'Illumina indices.')
    parser.add_argument('-n', '--num_indices', type=int, default=10)
    parser.add_argument('fastq_file', type=str)
    args = parser.parse_args()
    print_indices_summary(args.fastq_file, args.num_indices)

if __name__ == "__main__":
    main()