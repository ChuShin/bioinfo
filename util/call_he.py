#!/usr/bin/env python3

import argparse
import sys
import datetime
import numpy as np
import scipy.stats as stats
from collections import defaultdict
from itertools import chain

"""call_he.py: Given a coverage file, make HE calls."""



def load_covs(filename):
    """
    load coverage file
    :param filename:
    """

    covs = defaultdict(lambda: defaultdict(list))
    sample_metadata = {}

    linecount = 0
    with open(filename, 'r') as infile:
        for line in infile:
            dat = line.strip().split(',')
            linecount += 1
            if linecount <= 1:
                if dat[0] != "sample":
                    log_error(f"invalid file format: {line}")
            elif len(dat) == 7:
                covs[dat[0]][dat[1]].append([float(dat[4]), float(dat[6]), line])
    return covs

def normalize(covs):
    for sample in covs.keys():
        log_message(f"normalizing sample: {sample}")
        gene_A_covs = []
        gene_B_covs = []
        for chr_group in covs[sample].keys():
            for gene in covs[sample][chr_group]:
                gene_A_covs.append(gene[0])
                gene_B_covs.append(gene[1])
        observed_covs = [gene_A_covs, gene_B_covs]
        data = np.array(observed_covs)
        mean = np.mean(data)
        stdev = np.std(data)
        log_message(f" {sample} has mean= {mean:.2f}, stdev= {stdev:.2f}")

        z_scores = stats.zscore(data,axis=1)

        #one-sided filtering of z-score to remove extreme high cov.
        #filtered_entries = (z_scores < 3).all(axis=1)
        #disable for now.
        #new_df = df[filtered_entries]
        #print(f"{data}")
        #print(f"{z_scores}")

def log_error(message):
    """
    log error message with datetime and exit the program
    :param message:
    """

    dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{dt}]: ERROR {message}")
    sys.exit(1)

def log_message(message):
    """
    log message
    :param message:
    """

    dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{dt}]: {message}")


def main():
    parser = argparse.ArgumentParser(
        description="Given a coverage file, produce a summary of HE events.")
    parser.add_argument('-o', '--output_file', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('input_cov_file', type=str)
    args = parser.parse_args()
    covs = load_covs(args.input_cov_file)
    norm_cov = normalize(covs)

if __name__ == "__main__":
    main()
