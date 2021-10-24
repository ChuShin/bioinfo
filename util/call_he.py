#!/usr/bin/env python3

import argparse
import sys
from collections import namedtuple
from collections import defaultdict
from itertools import chain


"""call_he.py: Given a coverage file, make HE calls."""





def main():
    parser = argparse.ArgumentParser(
        description="Given a coverage file, produce a summary of HE events.")
    parser.add_argument('-o', '--output_file', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('input_cov_file', type=str)
    args = parser.parse_args()
    querys = load_query_assm_lengths(args.query_assm_info_file)
    covs = load_covs(args.input_cov_file)
#    summarize_simple(alns, querys,
#                     args.min_percent_id, args.min_aln_length)
#    summarize_combined(alns, querys, args.min_percent_id, args.min_aln_length)

    refs = load_ref_assm_lengths(args.query_assm_info_file)
    summarize_ref_simple(alns, refs,
                     args.min_percent_id, args.min_aln_length)

if __name__ == "__main__":
    main()
