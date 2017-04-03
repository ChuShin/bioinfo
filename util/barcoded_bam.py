#!/usr/bin/env python

import argparse
import pysam
import sys

def read_bam(filename, min_score):
    samfile = pysam.AlignmentFile(filename, 'rb')
    for read in samfile.fetch():
        read_tags = dict(read.tags)
        if 'MI' in read_tags and read_tags['AS'] >= min_score:
            print '%s\t%s\t%d\t%s\t%s' %(read.query_name, read.reference_name,
                                         read.reference_start, read_tags['MI'],
                                         read_tags['BX'])
    samfile.close()


def main():
    parser = argparse.ArgumentParser(
        description='Parser for barcoded BAM')
    parser.add_argument('-mi', '--molecule_only', type=bool, default=1)
    parser.add_argument('-ms', '--min_score', type=int, default=30)
    parser.add_argument('bam_filename', type=str)
    args = parser.parse_args()
    read_bam(args.bam_filename, args.min_score)


if __name__ == '__main__':
    main()