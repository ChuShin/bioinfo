#!/usr/bin/env python

import argparse
#import pysam
import sys
from collections import defaultdict

"""gean_choose_primary_models.py:
Given an input gene annotation file, choose the primary model (i.e. mRNA with
the longest ORF) and drop all alternate models."""




"""GFF class."""

class GFF:

    gff = defaultdict( lambda : defaultdict ( lambda : list ) )

    def __init__(self, line):
        self.line = line
        data = line.strip().split("\t")
        self.ref, self.src, self.type = data[:3]
        self.start = int(data[3])
        self.end = int(data[4])
        self.strand, self.score, self.phase, self.attributes = data[5:]



def read_gff(filename):
    with open(filename, 'r') as gff_file:
        for line in gff_file:
            gff = GFF(line)




def main():
    parser = argparse.ArgumentParser(
        description='Given an input gene annotation file, choose the primary '
                    'model (i.e. mRNA with the longest ORF) and drop all '
                    'alternate models. ')
    parser.add_argument('gff_filename', type=str)
    args = parser.parse_args()
    read_gff(args.gff_filename)


if __name__ == '__main__':
    main()



