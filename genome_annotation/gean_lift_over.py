#!/usr/bin/env python

import argparse
import pysam
import sys
from collections import defaultdict

"""gean_liftover.py:
Given an AGP file, convert coordinates in an input GFF/BED file between object
and component coordinates."""





def main():
    parser = argparse.ArgumentParser(
        description='Given an AGP file, convert coordinates in an input '
                    'GFF/BED file between object and component coordinates')
    parser.add_argument('gff_filename', type=str)
    parser.add_argument('bed_filename', type=str)
    parser.add_argument('to_component', type=bool, default=0)
    parser.add_argument('agp_filename', type=str)
    args = parser.parse_args()
    read_gff(args.gff_filename)

if __name__ == '__main__':
    main()
