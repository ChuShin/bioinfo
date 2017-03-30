#!/usr/bin/env python

import argparse
import sys
import re

def print_indices_summary(fastq_file, num_indices):
    with gzip.open(fastq_file, 'r') as infile:
        for line in infile:
            if line.startswith('@'):
                [header, index] = line.split(' ')
                print index



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