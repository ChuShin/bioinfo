#!/usr/bin/env python

import argparse
import pysam
import sys

def read_bam(filename):
    samfile = pysam.AlignmentFile(filename, 'rb')
    for read in samfile.fetch():
        read_tags = dict(read.tags)
        if read_tags['MI']:
            print '%s\t%s\t%d\t%s' %(read.query_name, read.reference_name,
                                     read.reference_start, read_tags['MI'])
    samfile.close()


def main():
    parser = argparse.ArgumentParser(
        description='Parser for barcoded BAM')
    parser.add_argument('-mi', '--molecule_only', type=boolean, default=TRUE)
    parser.add_argument('bam_filename', type=str)
    args = parser.parse_args()
    read_bam(args.bam_filename)


if __name__ == '__main__':
    main()