#!/usr/bin/python

import argparse
import re


"""calculate_n50.py:
Produce the assembly statistics given a WGA FASTA file."""


def read_sizes_from_fasta_file(filename, min_length):
    """Read FASTA sequence file return sort list and total length"""
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


def print_nties(n, cur, wrap, pad):
    """Print assembly stats"""
    print "%s value : %s %s %s size : %s" \
        %(n, cur, pad.rjust(wrap-len(cur)), n, "{:,}".format(sizes[j]))



def print_summary(fasta_file,min_length):
    [sizes, sum_length] = read_sizes_from_fasta_file(fasta_file, min_length)
    num_seq = "{:,}".format(len(sizes))
    max_len = "{:,}".format(sizes[0])
    min_len = "{:,}".format(sizes[-1])
    [count, f50, f60, f70, f80, f90] = [0, 0, 0, 0, 0, 0]
    wrap = 10
    pad = ' '

    print "Num sequences : %s %s Total bases : %s" \
        %(num_seq, pad.rjust(wrap-len(num_seq)), "{:,}".format(sum_length))

    print "Max. length   : %s %s Min. length : %s" \
        %(max_len, pad.rjust(wrap-len(max_len)), min_len)

    for j in range(0, len(sizes)):
        count+=sizes[j]
        cur = "{:,}".format(j+1)
        if count>=sum_length*0.5 and not(f50):
            print_nties("N50", cur, wrap, pad)
            f50 =1
        if count>=sum_length*0.6 and not(f60):
            print_nties("N60", cur, wrap, pad)
            f60 =1
        if count>=sum_length*0.7 and not(f70):
            print_nties("N70", cur, wrap, pad)
            f70 =1
        if count>=sum_length*0.8 and not(f80):
            print_nties("N80", cur, wrap, pad)
            f80 =1
        if count>=sum_length*0.9 and not(f90):
            print_nties("N90", cur, wrap, pad)
            f90 =1
            break


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