#!/usr/bin/env python

import argparse
import sys
from collections import namedtuple
from collections import defaultdict
from itertools import chain


"""anchor.py: Given an alignment file, produce a summary of where each query
sequence is landed."""


Alignment = namedtuple('Alignment', 'r_s, r_e, q_s, q_e, '
                                    'r_aln_len, q_aln_len, pid, r_len, '
                                    'q_len, r_cov, q_cov, ref, query')
flatten = chain.from_iterable


def merge_ranges(data):
    data = sorted(flatten(((start, 1), (stop, -1)) for start, stop in data))
    merged = []
    start, end, count = 0, 0, 0
    for pos, side in data:
        if count == 0:
            start = pos
        count += side
        if count == 0:
            end = pos
            merged.append([start, end])
    return merged


def load_coords(filename):
    alns = []
    with open(filename, 'r') as infile:
        for line in infile:
            dat = line.strip().split('\t')
            if len(dat) == 13:
                alns.append(Alignment._make(dat))
    return alns


def load_query_assm_lengths(filename):
    querys = defaultdict(lambda: defaultdict(int))
    with open(filename, 'r') as infile:
        for line in infile:
            key, val = line.strip().split('\t')
            query_name, c_name = key.split('_')
            querys[query_name][key] = float(val)
    return querys


def load_ref_assm_lengths(filename):
    refs = defaultdict(lambda: defaultdict(int))
    with open(filename, 'r') as infile:
        for line in infile:
            key, val = line.strip().split('\t')
            refs[key] = float(val)
    return refs


def summarize_ref_simple(alns, refs, min_pid, min_aln_len):
    """
    simple anchors considers only coverage in query, and makes no attempt to
    check for reference chaining
    :param alns:
    :param refs:
    :param min_pid:
    :param min_aln_len:
    """
    anchors = load_ref_anchors_simple(alns, min_pid, min_aln_len)
    for ref in anchors.keys():
        locs = anchors[ref]
#        print "\t", ref, contig, locs
        ranges = merge_ranges(locs)
#        print "\t\t", ranges
        total_bases = sum(e - s + 1 for s, e in ranges)
        print "%s\t%d\t%d\t%2.1f" \
                % (ref, total_bases, refs[ref],
                100 * total_bases / refs[ref]), ranges

def summarize_simple(alns, querys, min_pid, min_aln_len):
    """
    simple anchors considers only coverage in query, and makes no attempt to
    check for reference chaining
    :param alns:
    :param querys:
    :param min_pid:
    :param min_aln_len:
    """
    anchors = load_anchors_simple(alns, min_pid, min_aln_len)
    for query in anchors.keys():
        for ref in anchors[query].keys():
            for contig in anchors[query][ref].keys():
                locs = anchors[query][ref][contig]
#                print "\t", ref, contig, locs
                ranges = merge_ranges(locs)
#                print "\t\t", ranges
                total_bases = sum(e - s + 1 for s, e in ranges)
                print "%s\t%s\t%s\t%d/%d\t%2.1f" \
                      % (query, ref, contig, total_bases, querys[query][contig],
                         100 * total_bases / querys[query][contig]), locs



def summarize_combined(alns, querys, min_pid, min_aln_len):
    """
    simple anchors considers only coverage in query, and makes no attempt to
    check for reference chaining
    :param alns:
    :param querys:
    :param min_pid:
    :param min_aln_len:
    """
    anchors = load_anchors_combined(alns, min_pid, min_aln_len)
    for query in anchors.keys():
        for contig in anchors[query].keys():
            locs = anchors[query][contig]
#            print "\t", ref, contig, locs
            ranges = merge_ranges(locs)
#            print "\t\t", ranges
            total_bases = sum(e - s + 1 for s, e in ranges)
            print "%s\t%s\t%d/%d\t%2.1f" \
                % (query, contig, total_bases, querys[query][contig],
                    100 * total_bases / querys[query][contig]), locs


def load_anchors_combined(alns, min_pid, min_aln_len):
    """
    simple anchoring consider only coverage in query, no reference chaining
    :param alns:
    :param min_pid:
    :param min_aln_len:
    """
    anchors = defaultdict(lambda: defaultdict(list))
    for aln in alns:
        q_s, q_e = strandless(int(aln.q_s), int(aln.q_e))
        query_name = aln.query.split("_")[0]
        if float(aln.pid) >= min_pid and int(aln.q_aln_len) >= min_aln_len:
            anchors[query_name][aln.query].append([q_s, q_e])
    return anchors


def load_anchors_simple(alns, min_pid, min_aln_len):
    """
    simple anchoring consider only coverage in query, no reference chaining
    :param alns:
    :param min_pid:
    :param min_aln_len:
    """
    anchors = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for aln in alns:
        q_s, q_e = strandless(int(aln.q_s), int(aln.q_e))
        query_name = aln.query.split("_")[0]
        if float(aln.pid) >= min_pid and int(aln.q_aln_len) >= min_aln_len:
            anchors[query_name][aln.ref][aln.query].append([q_s, q_e])
    return anchors


def load_ref_anchors_simple(alns, min_pid, min_aln_len):
    """
    simple anchoring consider only coverage in query, no reference chaining
    :param alns:
    :param min_pid:
    :param min_aln_len:
    """
    anchors = defaultdict(list)
    for aln in alns:
        r_s, r_e = strandless(int(aln.r_s), int(aln.r_e))
        if float(aln.pid) >= min_pid and int(aln.r_aln_len) >= min_aln_len:
            anchors[aln.ref].append([r_s, r_e])
    return anchors

def strandless(aln_s, aln_e):
    if aln_s < aln_e:
        return aln_s, aln_e
    else:
        return aln_e, aln_s


def main():
    parser = argparse.ArgumentParser(
        description="Given an alignment file, produce a summary of where each "
                    "query sequence is landed.")
    parser.add_argument('-i', '--min_percent_id', type=int, default=99)
    parser.add_argument('-l', '--min_aln_length', type=int, default=1000)
    parser.add_argument('-b', '--query_assm_info_file', type=str)
    parser.add_argument('-o', '--output_file', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('input_coord_file', type=str)
    args = parser.parse_args()
    querys = load_query_assm_lengths(args.query_assm_info_file)
    alns = load_coords(args.input_coord_file)
#    summarize_simple(alns, querys,
#                     args.min_percent_id, args.min_aln_length)
#    summarize_combined(alns, querys, args.min_percent_id, args.min_aln_length)

    refs = load_ref_assm_lengths(args.query_assm_info_file)
    summarize_ref_simple(alns, refs,
                     args.min_percent_id, args.min_aln_length)

if __name__ == "__main__":
    main()
