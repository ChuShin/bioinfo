#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict

"""hic_assembly2agp.py:
Convert a Juicebox assembly file into AGP."""


# to-do: should create a separate class for agp file type
def read_agp_file(filename):
    chrs = defaultdict(lambda : defaultdict(lambda: list))
    with open(filename, 'r') as agp:
        for line in agp:
            if not line.startswith('#'):
                [object, object_beg, object_end, part_number,
                 component_type, component_id, component_beg,
                 component_end, strand] = line.strip().split('\t')
                if component_type == 'W' or component_type == 'D':
                    chrs[component_id][int(component_end)] = {'object': object,
                                            'object_beg': int(object_beg),
                                            'object_end': int(object_end),
                                            'component_beg' : int(component_beg),
                                            'component_end' :int(component_end),
                                            'strand' : strand}
    return chrs

def is_valid_feature(feature_name, feature_start, feature_end, component):
    if feature_start >= component['component_beg'] and \
        feature_start <= component['component_end']:
        if feature_end <= component['component_end']:
            return True
        else:
            print >> sys.stderr, 'Error: %s crossed component boundary' \
                                 %(feature_name)


def lookup(start, end, strand, component):
    # offset component start position
    new_start = start - component['component_beg']
    new_end = end - component['component_beg']
    new_strand = assign_strand(strand, component['strand'])
    if component['strand'] == '+' or component['strand'] == '?':
        new_start = new_start + component['object_beg']
        new_end = new_end + component['object_beg']
    elif component['strand'] == '-':
        tmp_pos = new_start
        new_start = component['object_end'] - new_end
        new_end = component['object_end'] - tmp_pos
    else:
        print >> sys.stderr, 'Error: unknown strand "%s" value, ' \
                             'please use (+,-,?)' %(component['strand'])
    return new_start, new_end, new_strand


def assign_strand(feature_strand, component_strand):
    if component_strand == '?':
        component_strand = '+'
    if feature_strand == component_strand:
        return '+'
    else:
        return '-'


def main():
    parser = argparse.ArgumentParser(
        description='Convert a Juicebox assembly file into AGP. ')
    parser.add_argument('-agp', '--agp_filename', type=str, default='default.agp')
    parser.add_argument('assm_filename', type=str)
    args = parser.parse_args()
    chrs = read_agp_file(args.agp_filename)
    if args.gff_filename is not None:
        gff_lift_over(args.gff_filename,chrs)
    elif args.bed_filename is not None:
        bed_lift_over(args.bed_filename,chrs)
    else:
        print >> sys.stderr, 'Please provide an input .assembly file'
        sys.exit(1)


if __name__ == '__main__':
    main()