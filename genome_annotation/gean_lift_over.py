#!/usr/bin/env python

import argparse
import pysam
import sys
from collections import defaultdict

"""gean_liftover.py:
Given an input GFF file ."""





def main():
    parser = argparse.ArgumentParser(
        description='Convert coordinates of GFF file onto new REFs '
                    'according to an AGP file. ')
    parser.add_argument('gff_filename', type=str)
    parser.add_argument('agp_filename', type=str)
    args = parser.parse_args()
    read_gff(args.gff_filename)

if __name__ == '__main__':
    main()
