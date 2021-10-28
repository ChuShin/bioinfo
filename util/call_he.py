#!/usr/bin/env python3

import argparse
import sys
import datetime
import numpy as np
import scipy.stats as stats
from collections import defaultdict
from itertools import chain
from itertools import islice

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
                covs[dat[0]][dat[1]].append([float(dat[4]), float(dat[6]), line.strip()])
    return covs

def normalize(covs):
    for sample in covs.keys():
        normalize_sample(sample, covs[sample])


def normalize_sample(sample, sample_cov):
    gene_A_covs = []
    gene_B_covs = []
    log_message(f"normalizing sample: {sample}")

    for chr_group in sample_cov.keys():
        for gene in sample_cov[chr_group]:
            gene_A_covs.append(gene[0])
            gene_B_covs.append(gene[1])
    observed_covs = [gene_A_covs, gene_B_covs]
    data = np.array(observed_covs)
    mean = np.mean(data)
    stdev = np.std(data)
    log_message(f"sample= {sample}, mean= {mean:.2f}, stdev= {stdev:.2f}")
    z_scores = stats.zscore(data,axis=1)
    call_he_events(z_scores, sample_cov)

    #one-sided filtering of z-score to remove extreme high cov.
    #filtered_entries = (z_scores < 3).all(axis=1)
    #disable for now.
    #new_df = df[filtered_entries]
    #print(f"{data}")
    #print(f"{z_scores}")


def call_he_events(z_scores, sample_cov):
    gene_count = len(z_scores[0])
    he_array = []
    for i in range(gene_count):
        zscore_A = z_scores[0][i]
        zscore_B = z_scores[1][i]
        he_type = call_he_type(zscore_A, zscore_B)
        he_array.append(he_type)
    seeds = find_seeds(he_array)
    he_events = extend_seeds(seeds)
    print_he_output(he_events, z_scores, sample_cov)


def print_he_output(he_events, z_scores, sample_cov):
    for chr_group in sample_cov.keys():
        for i, gene in enumerate(sample_cov[chr_group]):
            print(f"{gene[2]},{z_scores[0][i]},{z_scores[1][i]},{he_events[i]}")




def find_seeds(he_array):
    seed_size = 3
    gene_count = len(he_array) - seed_size
    he_seeds = [''] * len(he_array)
    for i, he_type in enumerate(he_array):
        #print(f"{i} {he_type}")
        if i <= gene_count:
            if is_hit(he_array, i, seed_size):
                #print(f"{i} is hit")
                for j in range(i,i+seed_size):
                    he_seeds[j] = he_type
#    print(f"{he_events}")
    return he_seeds

def extend_seeds(he_array):
    #skip
    return he_array


def is_hit(he_array, i, seed_size):
    if he_array[i] == "norm/norm":
        return False
    for j in range(i, i+seed_size):
        #print(f"j is {j}")
        if he_array[i] != he_array[j]:
            return False
    return True


def call_he_type(zscore_A, zscore_B):
    cutoff = 1.5
    if zscore_A < -1*cutoff and zscore_B < -1*cutoff:
        return "del/del"
    if zscore_A < -1*cutoff and abs(zscore_B) < cutoff:
        return "del/norm"
    if zscore_A < -1*cutoff and zscore_B > cutoff:
        return "del/dup"
    if abs(zscore_A) < cutoff and zscore_B < -1*cutoff:
        return "norm/del"
    if abs(zscore_A) < cutoff and abs(zscore_B) < cutoff:
        return "norm/norm"
    if abs(zscore_A) < cutoff and zscore_B > cutoff:
        return "norm/dup"
    if zscore_A > cutoff and zscore_B < -1*cutoff:
        return "dup/del"
    if zscore_A > cutoff and abs(zscore_B) < cutoff:
        return "dup/norm"
    if zscore_A > cutoff and zscore_B > cutoff:
        return "dup/dup"
    return "undefined"


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
