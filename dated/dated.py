#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import asyncio
import tempfile
import shutil
from multiprocessing import Process, Pool
from functools import partial
from Bio import SeqIO
from Bio.Phylo.PAML import codeml
from Bio.Align.Applications import ClustalwCommandline

"""DATED : DATing Evolutionary events and Divergence times."""


def read_sequence(filename):
    return SeqIO.to_dict(SeqIO.parse(filename, "fasta"))

def read_pairs(filename):
    ortho_pairs = []
    with open(filename, 'r') as pairs:
        for pair in pairs:
            ortho_pairs.append(pair.rstrip().split(","))
    return ortho_pairs

def calculate_ds(pairs, pep_seqdb, cds_seqdb):
    ds_func = partial(get_ds, pep_seqdb, cds_seqdb)
    with Pool(processes=8) as pool:
        pool.map(ds_func, pairs)


def get_ds(pep_seqdb, cds_seqdb, pairs):
#    for [seqA, seqB] in pairs:
        [seqA, seqB] = pairs
        try:
            tmp_folder = tempfile.mkdtemp(dir="./tmp")
            cwd = os.getcwd()
            os.chdir(tmp_folder)
            pep_seq_path = "pep_pair.fasta"
            cds_seq_path = "cds_pair.fasta"
            aln_path = "pep_pair.aln"
            pal_path = "pep_pair.pal2nal"
            run_clustalw(seqA, seqB, pep_seqdb, pep_seq_path)
            run_pal2nal(seqA, seqB, cds_seqdb, cds_seq_path,
                        aln_path, pal_path)
            ds = run_codeml(tmp_folder, pal_path)
            print(",".join(map(str,(seqA, seqB, ds))))
        finally:
            try:
                os.chdir(cwd)
                shutil.rmtree(tmp_folder)  # delete directory
            except OSError as exc:
                raise exc



def run_clustalw(seqA, seqB, seqdb, seq_path):
    with open(seq_path, "w") as seqfile:
        SeqIO.write(seqdb[seqA], seqfile, "fasta")
        SeqIO.write(seqdb[seqB], seqfile, "fasta")
    clustalw = ClustalwCommandline('clustalw', infile=seq_path)
    stdout, stderr = clustalw()

def run_pal2nal(seqA, seqB, seqdb, cds_seq_path, aln_path, pal_path):
    with open(cds_seq_path, "w") as seqfile:
        SeqIO.write(seqdb[seqA], seqfile, "fasta")
        SeqIO.write(seqdb[seqB], seqfile, "fasta")
    f = open(pal_path, "w")
    subprocess.run(["pal2nal.pl",aln_path, cds_seq_path,
                    "-nogap", "-output", "paml"], stdout=f)

def run_codeml(tmp_folder, pal_path):
    cml = codeml.Codeml(alignment=pal_path,
                        out_file="pair.ks",
                        tree="pep_pair.dnd")
    cml.read_ctl_file("../../config/codeml.ctl")
    results = cml.run().get("pairwise")
    prot1 = next(iter(results.values()))
    for prot2, attributes in prot1.items():
        ds = (attributes.get("dS"))
    return ds


def main():
    parser = argparse.ArgumentParser(
        description='Given a list of ')
    parser.add_argument('pep_filename', type=str)
    parser.add_argument('cds_filename', type=str)
    parser.add_argument('pair_list', type=str)
    args = parser.parse_args()
    pep_seqdb = read_sequence(args.pep_filename)
    cds_seqdb = read_sequence(args.cds_filename)
    pairs = read_pairs(args.pair_list)
    calculate_ds(pairs, pep_seqdb, cds_seqdb)




if __name__ == '__main__':
    main()