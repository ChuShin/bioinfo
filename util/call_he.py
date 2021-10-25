#!/usr/bin/env python3

import argparse
import sys
import datetime
from collections import namedtuple
from collections import defaultdict
from numpy import genfromtxt
from itertools import chain
from datetime import date

"""call_he.py: Given a coverage file, make HE calls."""

Coverage = namedtuple('Coverage', 'sample, group, pos, geneA, covA, geneB, covB')
flatten = chain.from_iterable


def load_covs(filename):
    covs = []
    with open(filename, 'r') as infile:
        for line in infile:
            dat = line.strip().split(',')
            if len(covs) == 0:
                if dat[0] != "sample":
                    log_error("invalid file format.")
            if len(dat) == 7:
                covs.append(Coverage._make(dat))
    return covs


def log_error(message):
    dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{dt}]: {message}")
    sys.exit(1)

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
