#!/usr/bin/python

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--min_length", type=int, default=200, 
        help="only consider sequences greater than this length")
parser.add_argument("fasta_file", type=str)
args = parser.parse_args()

def print_nties(n):
    """Print assembly stats"""
    print "%s value : %s %s %s size : %s" \
        %(n, cur, pad.rjust(wrap-len(cur)), n, "{:,}".format(sizes[j]))

def read_sizes_from_fasta_file(filename):
    """Read FASTA sequence file return sort list and total lenght"""
    f = open(filename, 'r')
    line = f.readline()
    sizes = []
    sum_length = 0
    cur_length = 0

    while line:
        line=line.rstrip("\n")
        if line.startswith(">"):
            if cur_length > args.min_length:
                sizes.append(cur_length)
                sum_length += cur_length
            cur_length = 0
        else:
            cur_length += len(line)

        line = f.readline()

    if cur_length > args.min_length:
        sizes.append(cur_length)
        sum_length += cur_length

    sizes.sort(reverse=True)

    return sizes, sum_length

[sizes, sum_length] = read_sizes_from_fasta_file(args.fasta_file)
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
        print_nties("N50")
        f50 =1
    if count>=sum_length*0.6 and not(f60):
        print_nties("N60")
        f60 =1
    if count>=sum_length*0.7 and not(f70):
        print_nties("N70")
        f70 =1
    if count>=sum_length*0.8 and not(f80):
        print_nties("N80")
        f80 =1
    if count>=sum_length*0.9 and not(f90):
        print_nties("N90")
        f90 =1
        break
