#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
import argparse
import sys
import pysam

HI_CONF = 95.0

def load_pos_in_bed_file(pos_file):
    queries = {}
    with open(pos_file, 'r') as infile:
        for line in infile:
            q = line.strip().split('\t')
            queries[q[0]] = [q[1], q[2]]
    return queries

def get_confidence_level(perc):
    if perc >= HI_CONF:
        return '100'
    elif perc < HI_CONF:
        return '50'

def parse_sam(input_sam_file, output_file, queries):
    bam_file = pysam.Samfile(input_sam_file, "rb");
    gff_file = open(output_file, 'w')
    log_file = open(output_file+'.log', 'w')
    for read in bam_file:
        if not read.is_unmapped:   #if it's mapped
            pos_on_query, alleles = queries[read.query_name]
            pos_on_query = int(pos_on_query) - 1    # pysam is zero-based
            ref_positions = read.get_reference_positions(full_length=True)
            strand = '+'
            if read.is_reverse:
                pos_on_query = read.query_length - pos_on_query - 1
                strand = '-'
            ref_snp_pos = ref_positions[pos_on_query]
            prop_aligned =  (( len(read.get_reference_positions()) -
                read.get_tag('NM'))  / read.query_length ) * 100
            if ref_snp_pos != None:
                ## gff output 1-based
                ref_snp_pos = ref_snp_pos + 1
                gff_line = "%s\tUofS\tmarker\t%d\t%d\t%s\t%s\t.\t%s\n" %(
                    bam_file.getrname(read.reference_id), ref_snp_pos,
                    ref_snp_pos, get_confidence_level(prop_aligned), strand,
                    'Name='+read.query_name+';Note=PID:'+'%.2f' % prop_aligned)
                gff_file.write(gff_line)
            elif ref_snp_pos == None:
                log_line = "[No SNP position error]\t%s\t%s\t%s\t%s\t%s\n" %(
                    read.query_name, pos_on_query, alleles,
                    bam_file.getrname(read.reference_id), read.cigarstring)
                log_file.write(log_line)
                ## log_file.write(str(ref_positions)+'\n')   ## debug

def main():
    parser = argparse.ArgumentParser(
        description="Given a read alignment file and a base position on reads, "
            "produce a GFF file of where each base pos is located on the reference.")
    parser.add_argument('-b', '--base_pos_info_file', type=str)
    parser.add_argument('-o', '--output_file', type=str)
    parser.add_argument('input_sam_file', type=str)
    args = parser.parse_args()
    queries = load_pos_in_bed_file(args.base_pos_info_file)
    parse_sam(args.input_sam_file, args.output_file, queries)


if __name__ == "__main__":
    main()

