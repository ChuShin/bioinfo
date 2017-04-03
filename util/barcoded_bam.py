#!/usr/bin/env python

import argparse
import pysam
import sys

def read_bam(filename):
    samfile = pysam.AlignmentFile(in_file, 'rb')
    for read in samfile.fetch():
     if read.is_proper_pair:
             pairedreads.write(read)

pairedreads.close()
samfile.close()


def main():
    parser = argparse.ArgumentParser(
        description='Parser for barcoded BAM')
    parser.add_argument('-mi', '--molecule_only', type=boolean, default=TRUE)
    parser.add_argument('bam_filename', type=str)


if __name__ == '__main__':
    main()