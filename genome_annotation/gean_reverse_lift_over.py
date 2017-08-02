#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict

"""gean_reverse_lift_over.py:
Given an AGP file, convert coordinates in an input GFF/BED file from object
to component coordinates."""


# to-do: should create a separate class for agp file type
def read_agp_file(filename):
    chrs = defaultdict(lambda : defaultdict(lambda: list))
    with open(filename, 'r') as agp:
        for line in agp:
            if not line.startswith('#'):
                [object, object_beg, object_end, part_number,
                 component_type, component_id, component_beg,
                 component_end, strand] = line.strip().split('\t')
                if component_type == 'W':
                    chrs[component_id] = {'object': object,
                                          'object_beg': int(object_beg),
                                          'object_end': int(object_end),
                                          'strand' : strand}
    return chrs

# to-do: should create a separate class for bed file type
def lift_over(filename, chrs):
    with open(filename, 'r') as bed:
        for line in bed:
            try:
                [chr, start, end, feature, score, strand] = \
                    line.strip().split('\t')
                component = chrs[chr]
                new_start, new_end, new_strand = \
                    lookup(int(start), int(end), strand, component)
                print '%s\t%d\t%d\t%s\t%s\t%s' %(
                    component['object'], new_start, new_end,
                    feature, score, new_strand)
            except ValueError:
                print >> sys.stderr, 'Invalid 6-col bed file format.'
                sys.exit(1)
            except Exception, e:
                print >> sys.stderr, 'Exception: %s' %str(e)
                sys.exit(1)

def lookup(start, end, strand, component):
    new_strand = assign_strand(strand, component['strand'])
    if component['strand'] == '+' or component['strand'] == '?':
        new_start = component['object_beg'] + start - 1
        new_end = component['object_beg'] + end - 1
    elif component['strand'] == '-':
        new_start = component['object_end'] - end + 1
        new_end = component['object_end'] - start + 1
    return new_start, new_end, new_strand


def assign_strand(feature_strand, component_strand):
    if component_strand == '?':
        component_strand = '+'
    if feature_strand == component_strand:
        return '+'
    elif feature_strand != component_strand:
        return '-'



def main():
    parser = argparse.ArgumentParser(
        description='Given an AGP file, convert coordinates in an input '
                    'GFF/BED file from object to component coordinates')
#    parser.add_argument('-gff', '--gff_filename', type=str)
    parser.add_argument('-bed', '--bed_filename', type=str)
    parser.add_argument('-c', '--to_component', type=bool, default=0)
    parser.add_argument('agp_filename', type=str)
    args = parser.parse_args()
    chrs = read_agp_file(args.agp_filename)
    lift_over(args.bed_filename,chrs)

if __name__ == '__main__':
    main()