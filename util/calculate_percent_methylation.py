#!/usr/bin/env python

import argparse
import sys
from collections import namedtuple
from collections import defaultdict

"""calculate_percent_methylation.py:
Given an intersected BED file produced by bedtools-intersectbed,
calculate the methylation percentage for each feature (column4)."""


Methyl_Feature = namedtuple('methyl_feature', 'chr, start, end, feature, score, '
                                    'strand, bs_chr, bs_start, bs_end, '
                                    'bs_context, bs_type, bs_strand, overlap')


def load_intersections(intersect_file):
    methyl = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    with open(intersect_file, 'r') as infile:
        for line in infile:
            data = line.strip().split('\t')
            coord = '\t'.join(data[0:4])
            intersection = Methyl_Feature._make(data)
            if int(intersection.overlap) > 0:
                methyl[coord][intersection.bs_context][intersection.bs_type]+=1
    return methyl

def calculate_percent_methylation(intersections, min_obs_nuc):
    for coord in intersections:
        for context in intersections[coord]:
            mc = intersections[coord][context]['5mC']
            uc = intersections[coord][context]['C']
            if (uc+mc) > min_obs_nuc:
                print '%s\t%s\t%d\t%d\t%.1f' %(coord, context, mc, uc+mc,
                                               mc/(mc+uc)*100)





def main():
    parser = argparse.ArgumentParser(
        description='Given an intersect bed output file, produce a summary of '
                    '%methyl for each input feature.')
    parser.add_argument('-m', '--min_obs_nuc', type=int, default=1)
    parser.add_argument('-o', '--output_file', type=argparse.FileType('w'),
                        default=sys.stdout)
    #parser.add_argument('intersect_file', type=str)
    args = parser.parse_args()
    intersections = load_intersections('../test_data/methylation_bed.txt')
    calculate_percent_methylation(intersections, args.min_obs_nuc)

if __name__ == "__main__":
    main()