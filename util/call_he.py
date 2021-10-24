#!/usr/bin/env python3

import argparse
import sys
from collections import namedtuple
from collections import defaultdict
from numpy import genfromtxt
from itertools import chain


"""call_he.py: Given a coverage file, make HE calls."""


def load_covs(filename):
    covs = []
#    with open(filename, 'r') as infile:
#        for line in infile:
#            dat = line.strip().split(',')
#            if len(dat) == 13:
#                alns.append(Alignment._make(dat))
    covs = genfromtxt(filename, delimiter=',')
    print covs
    return covs



def main():
    parser = argparse.ArgumentParser(
        description="Given a coverage file, produce a summary of HE events.")
    parser.add_argument('-o', '--output_file', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('input_cov_file', type=str)
    args = parser.parse_args()
    covs = load_covs(args.input_cov_file)

if __name__ == "__main__":
    main()
