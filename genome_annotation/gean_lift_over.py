#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict

"""gean_lift_over.py:
Given an AGP file, convert coordinates in an input GFF/BED file from component
coordinates into object coordinates."""


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
                    chrs[component_id][int(component_end)] = {'object': object,
                                            'object_beg': int(object_beg),
                                            'object_end': int(object_end),
                                            'component_beg' : int(component_beg),
                                            'component_end' :int(component_end),
                                            'strand' : strand}
    return chrs

# to-do: should create a separate class for bed file type
def bed_lift_over(filename, chrs):
    with open(filename, 'r') as bed:
        for line in bed:
            try:
                [chr, start, end, feature, score, strand] = \
                    line.strip().split('\t')
                if chr in chrs.keys():
                    for position in sorted(chrs[chr]):
                        if position >= int(start):
                            component = chrs[chr][position]
                            new_start, new_end, new_strand = \
                                lookup(int(start), int(end), strand,
                                        component)
                            new_start, new_end, new_strand = \
                                lookup(int(start), int(end), strand, component)
                            print '%s\t%d\t%d\t%s\t%s\t%s' %(
                                component['object'], new_start, new_end,
                                feature, score, new_strand)
                            break
                else:
                    print '%s does not exists in AGP' %(chr)
            except ValueError:
                print >> sys.stderr, 'Invalid 6-col bed file format.'
                sys.exit(1)
            except Exception, e:
                print >> sys.stderr, 'Exception: %s' %str(e)
                sys.exit(1)


# to-do: should create a separate class for bed file type
def gff_lift_over(filename, chrs):
    with open(filename, 'r') as gff:
        for line in gff:
            try:
                if line.startswith('#'):
                    print line
                else:
                    [chr, source, feature_type, start, end, score, strand,
                     frame, feature_name] = line.strip().split('\t')
                    if chr in chrs.keys():
                        for position in sorted(chrs[chr]):
                            component = chrs[chr][position]
                            if is_valid_feature(feature_name, int(start),
                                                int(end), component):
                                new_start, new_end, new_strand = \
                                    lookup(int(start), int(end), strand,
                                           component)
                                print '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s' %(
                                    component['object'], source, feature_type,
                                    new_start, new_end, score, strand, frame,
                                    feature_name)
                                break
                    else:
                        print >> sys.stderr, '%s does not exists in AGP' %(chr)
                        sys.exit(1)
            except ValueError:
                print >> sys.stderr, 'Invalid 6-col bed file format.'
                sys.exit(1)
            except Exception, e:
                print >> sys.stderr, 'Exception: %s' %str(e)
                sys.exit(1)

def is_valid_feature(feature_name, feature_start, feature_end, component):
    if feature_start >= component['component_end']:
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
        description='Given an AGP file, convert coordinates in an input '
                    'GFF/BED file from object to component coordinates')
    parser.add_argument('-gff', '--gff_filename', type=str)
    parser.add_argument('-bed', '--bed_filename', type=str)
    parser.add_argument('-c', '--to_component', type=bool, default=0)
    parser.add_argument('agp_filename', type=str)
    args = parser.parse_args()
    chrs = read_agp_file(args.agp_filename)
    if args.gff_filename is not None:
        gff_lift_over(args.gff_filename,chrs)
    elif args.bed_filename is not None:
        bed_lift_over(args.bed_filename,chrs)
    else:
        print >> sys.stderr, 'Please provide either a GFF or BED file'
        sys.exit(1)


if __name__ == '__main__':
    main()