#!/usr/bin/python

import argparse
import re


"""calculate_n50.py:
Produce the assembly statistics given a WGA FASTA file."""


def read_sizes_from_fasta_file(filename, min_length):
    """Return sort sizes and total WGA length (>min_length) in FASTA file"""
    f = open(filename, 'r')
    line = f.readline()
    sizes = []
    sum_length = 0
    cur_length = 0

    while line:
        line=line.rstrip("\n")
        if line.startswith(">"):
            if cur_length >= min_length:
                sizes.append(cur_length)
                sum_length += cur_length
            cur_length = 0
        else:
            cur_length += len(line)

        line = f.readline()

    if cur_length >= min_length:
        sizes.append(cur_length)
        sum_length += cur_length

    sizes.sort(reverse=True)

    return sizes, sum_length


def print_nties(nties, cur, size, wrap, pad):
    """Print assembly stats"""
    print "N%d value : %s %s %s size : %s" \
        %(nties, cur, pad.rjust(wrap-len(cur)), nties, "{:,}".format(size))


def print_summary(fasta_file,min_length):
    [sizes, sum_length] = read_sizes_from_fasta_file(fasta_file, min_length)
    num_seq = "{:,}".format(len(sizes))
    max_len = "{:,}".format(sizes[0])
    min_len = "{:,}".format(sizes[-1])
    count = 0
    wrap = 10
    pad = ' '

    print "Num sequences : %s %s Total bases : %s" \
        %(num_seq, pad.rjust(wrap-len(num_seq)), "{:,}".format(sum_length))

    print "Max. length   : %s %s Min. length : %s" \
        %(max_len, pad.rjust(wrap-len(max_len)), min_len)

    for j in range(0, len(sizes)):
        count+=sizes[j]
        cur = "{:,}".format(j+1)
        for nties in range(50, 90, 10):
            if count >= sum_length*nties/100:
                print_nties(nties, cur, sizes[j], wrap, pad, )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_length", type=int, default=200,
        help="only consider sequences greater than this length")
    parser.add_argument('-l', '--min_seq_length', type=int, default=200)
#    parser.add_argument("fasta_file", type=str)
    args = parser.parse_args()
    print_summary('../test_data/wga.fasta',args.min_seq_length)


if __name__ == "__main__":
    main()