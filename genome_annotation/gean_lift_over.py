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
            if line.startswith('#'):
                pass
            [object, object_beg, object_end, part_number,
             component_type, component_id, component_beg,
             component_end, strand] = line.strip().split('\t')
            if component_type == 'W':
                chrs[object][object_end] = {'object_beg': int(object_beg),
                                            'object_end': int(object_end),
                                            'component_id': component_id,
                                            'component_beg' : int(component_beg),
                                            'component_end' :int(component_end),
                                            'strand' : strand}
    return chrs

# to-do: should create a separate class for bed file type
def lift_over(filename, chrs):
    with open(filename, 'r') as bed:
        for line in bed:
            try:
                [chr, start, end, feature] = \
                    line.strip().split('\t')
                for position in chrs[chr]:
                    if position >= start:
                        component = chrs[chr][position]
                        new_start, new_end = lookup(int(start), int(end),
                                                    component)
                        print '%s\t%d\t%d\t%s' %(component['component_id'],
                                                 new_start, new_end, feature)
                        break
            except ValueError:
                print >> sys.stderr, 'Invalid 4-col bed file, please check format.'
                sys.exit(1)
            except Exception, e:
                print >> sys.stderr, 'Exception: %s' %str(e)
                sys.exit(1)

def lookup(start, end, component):
    if component['object_end'] < end:
        print >> sys.stderr, 'Feature [%d, %d] crossed an object boundary %s' \
                             %(start, end, component)
    new_start = start - component['object_beg']
    new_end = end - component['object_beg']
    if component['strand'] == '+':
        new_start = new_start + component['component_beg']
        new_end = new_end + component['component_beg']
    elif component['strand'] == '-':
        tmp_pos = new_start
        new_start = component['object_end'] - new_end
        new_end = component['object_end'] - tmp_pos
    return new_start, new_end

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
