#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict

"""gean_liftover.py:
Given an AGP file, convert coordinates in an input GFF/BED file between object
and component coordinates."""


# to-do: should create a separate class for agp file type
def read_agp_file(filename):
    chrs = defaultdict(lambda : defaultdict(lambda: list))
    with open(filename, 'r') as agp:
        for line in agp:
            [object, object_beg, object_end, part_number,
             component_type, component_id, component_beg,
             component_end, strand] = line.strip().split('\t')
            if component_type == 'W':
                chrs[object][object_beg] = [component_id, component_beg,
                                            component_end, strand]
    return chrs

# to-do: should create a separate class for bed file type
def lift_over(filename, chrs):
    with open(filename, 'r') as bed:
        for line in bed:
            try:
                [chr, start, end, feature, score, strand] = \
                    line.strip().split('\t')


            except ValueError:
                print >> sys.stderr, 'Invalid 6-col bed file, please check format.'
                sys.exit(1)
            except Exception, e:
                print >> sys.stderr, 'Exception: %s' %str(e)
                sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Given an AGP file, convert coordinates in an input '
                    'GFF/BED file between object and component coordinates')
#    parser.add_argument('-gff', '--gff_filename', type=str)
    parser.add_argument('-bed', '--bed_filename', type=str)
    parser.add_argument('-c', '--to_component', type=bool, default=0)
    parser.add_argument('agp_filename', type=str)
    args = parser.parse_args()
    chrs = read_agp_file(args.agp_filename)
    lift_over(args.bed_filename,chrs)

if __name__ == '__main__':
    main()
