import argparse
from collections import defaultdict
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from standard_logger import get_logger
from tandem_aligner_prototype import get_rare_kmers
from utils.os_utils import expandpath, smart_makedirs
from utils.git import get_git_revision_short_hash
from utils.bio import read_bio_seqs


SCRIPT_FN = os.path.realpath(__file__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--asm1", help="Asm1", required=True)
    parser.add_argument("--asm2", help="Asm2", required=True)
    parser.add_argument("--asm1-name", default="")
    parser.add_argument("--asm2-name", default="")
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("-k", type=int, default=7)
    parser.add_argument("--rare", type=int, default=1)
    params = parser.parse_args()

    params.asm1 = expandpath(params.asm1)
    params.asm2 = expandpath(params.asm2)
    if not os.path.isfile(params.asm1):
        logger.error(f'File does not exists --asm1 == {params.asm1}')
        sys.exit(1)
    if not os.path.isfile(params.asm2):
        logger.error(f'File does not exists --asm2 == {params.asm2}')
        sys.exit(1)

    params.outdir = expandpath(params.outdir)
    smart_makedirs(params.outdir)
    return params


def get_kmer2coords(s, kmers):
    k = len(next(iter(kmers)))
    kmer2coords = defaultdict(list)
    for i in range(len(s)-k+1):
        kmer = s[i:i+k]
        if kmer in kmers:
            kmer2coords[kmer].append(i)
    return kmer2coords


def get_pairs(kmers2coords1, kmers2coords2):
    pairs = []
    for kmer, coords1 in kmers2coords1.items():
        k = len(kmer)
        for coord2 in kmers2coords2[kmer]:
            for coord1 in coords1:
                pairs.append((coord1, coord2))
                pairs.append((coord1 + k, coord2 + k))
                pairs.append((np.nan, np.nan))
    return pairs


def main():
    params = parse_args()

    logfn = os.path.join(params.outdir, 'dotplot.log')
    global logger
    logger = get_logger(logfn,
                        logger_name='dotplot')
    logger.info(f'Constructing dotplot with {SCRIPT_FN} started')
    logger.info('cmd: {}'.format(sys.argv))
    logger.info('git hash: {}'.format(get_git_revision_short_hash()))

    asm1 = read_bio_seqs(params.asm1)
    asm2 = read_bio_seqs(params.asm2)

    asm1_name, asm1 = next(iter(asm1.items()))
    asm2_name, asm2 = next(iter(asm2.items()))
    asm1_name, asm2_name = asm1_name.split('|')[0], asm2_name.split('|')[0]

    if params.asm1_name != "":
        asm1_name = params.asm1_name
    if params.asm2_name != "":
        asm2_name = params.asm2_name


    rare_kmers1 = get_rare_kmers(asm1, params.k, params.rare)
    rare_kmers2 = get_rare_kmers(asm2, params.k, params.rare)
    rare_kmers = rare_kmers1 & rare_kmers2

    kmer2coords1 = get_kmer2coords(asm1, rare_kmers)
    kmer2coords2 = get_kmer2coords(asm2, rare_kmers)

    pairs = get_pairs(kmer2coords1, kmer2coords2)
    plt.figure(figsize=(5, 5), dpi=200)
    plt.plot([pair[0] for pair in pairs],
             [pair[1] for pair in pairs],
             marker=".", ms=0.5, alpha=1, color='black')
    plt.xlabel(asm1_name)
    plt.ylabel(asm2_name)
    plt.xlim(0)
    plt.ylim(0)
    plt.grid(b=True, which='major', axis='both')
    plt.xticks(rotation=30)
    plt.title(f'k={params.k}, MaxRareCount={params.rare}')
    outfn = os.path.join(params.outdir, 'dotplot.png')
    plt.savefig(outfn, format='png', bbox_inches="tight")
    plt.close()



if __name__ == "__main__":
    main()
