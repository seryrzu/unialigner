import argparse
from collections import Counter
import os
import sys


import numpy as np

from standard_logger import get_logger
from utils.os_utils import expandpath, smart_makedirs
from utils.git import get_git_revision_short_hash
from utils.bio import read_bio_seqs


SCRIPT_FN = os.path.realpath(__file__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--asm1", help="Asm1", required=True)
    parser.add_argument("--asm2", help="Asm2", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--k", type=int, default=7)
    parser.add_argument("--rare", type=int, default=1)
    params = parser.parse_args()

    params.asm1 = expandpath(params.asm1)
    params.asm2 = expandpath(params.asm2)
    if not os.path.isfile(params.asm1):
        logger.error(f"File does not exists --asm1 == {params.asm1}")
        sys.exit(1)
    if not os.path.isfile(params.asm2):
        logger.error(f"File does not exists --asm2 == {params.asm2}")
        sys.exit(1)

    params.outdir = expandpath(params.outdir)
    smart_makedirs(params.outdir)
    return params


def get_rare_kmers(s, k, rare):
    kmers_cnt = Counter([s[i : i + k] for i in range(len(s) - k + 1)])
    rare_kmers = set(kmer for kmer, cnt in kmers_cnt.items() if cnt <= rare)
    return rare_kmers


def get_coords(s, kmers):
    k = len(next(iter(kmers)))
    s_tr = []
    for i in range(len(s) - k + 1):
        kmer = s[i : i + k]
        if kmer in kmers:
            s_tr.append(i)
    return s_tr


def glob_align(s1, s2, coords_1, coords_2, k, match=1, mismatch=-np.inf, indel=0):
    n, m = len(coords_1) + 1, len(coords_2) + 1
    score = [[0 for j in range(m)] for i in range(n)]

    for i in range(1, n):
        score[i][0] = score[i - 1][0] + indel
    for j in range(1, m):
        score[0][j] = score[0][j - 1] + indel

    for i in range(1, n):
        for j in range(1, m):
            x = coords_1[i - 1]
            y = coords_2[j - 1]

            kmer1 = s1[x : x + k]
            kmer2 = s2[y : y + k]
            add = match if kmer1 == kmer2 else mismatch
            score[i][j] = max(
                score[i - 1][j - 1] + add,
                score[i - 1][j] + indel,
                score[i][j - 1] + indel,
            )
    a1, a2 = [], []
    i, j = n - 1, m - 1
    while i > 0 and j > 0:
        x = coords_1[i - 1]
        y = coords_2[j - 1]

        kmer1 = s1[x : x + k]
        kmer2 = s2[y : y + k]
        add = match if kmer1 == kmer2 else mismatch
        if score[i][j] == score[i - 1][j - 1] + add:
            a1.append((i - 1, x))
            a2.append((j - 1, y))
            i -= 1
            j -= 1
        elif score[i][j] == score[i - 1][j] + indel:
            a1.append((i - 1, x))
            a2.append("-")
            i -= 1
        else:
            assert score[i][j] == score[i][j - 1] + indel
            a1.append("-")
            a2.append((j - 1, y))
            j -= 1

    a1 += ["-"] * j
    while j > 0:
        a2.append((j - 1, coords_2[j - 1]))
        j -= 1

    a2 += ["-"] * i
    while i > 0:
        a1.append((i - 1, coords_2[i - 1]))
        i -= 1

    a1 = a1[::-1]
    a2 = a2[::-1]

    return a1, a2


def align(asm1_fn, asm2_fn, outdir, n_threads, k, rare):
    asm1 = read_bio_seqs(asm1_fn)
    asm2 = read_bio_seqs(asm2_fn)
    assert len(asm1) == 1 and len(asm2) == 1
    asm1_name, asm1 = list(asm1.items())[0]
    asm2_name, asm2 = list(asm2.items())[0]
    rare_kmers_1 = get_rare_kmers(asm1, k, rare)
    rare_kmers_2 = get_rare_kmers(asm2, k, rare)
    rare_kmers = rare_kmers_1 & rare_kmers_2
    coords_1 = get_coords(asm1, rare_kmers)
    coords_2 = get_coords(asm2, rare_kmers)
    a1, a2 = glob_align(asm1, asm2, coords_1, coords_2, k)

    align_fn = os.path.join(outdir, "glob_align.txt")
    with open(align_fn, "w") as f:
        for c1, c2 in zip(a1, a2):
            print(c1, c2, file=f)


def main():
    params = parse_args()
    logfn = os.path.join(params.outdir, "align_centromeres.log")
    global logger
    logger = get_logger(logfn, logger_name="tandem_aligner")
    logger.info(f"Aligning with {SCRIPT_FN} started")
    logger.info("cmd: {}".format(sys.argv))
    logger.info("git hash: {}".format(get_git_revision_short_hash()))

    align(
        asm1_fn=params.asm1,
        asm2_fn=params.asm2,
        outdir=params.outdir,
        n_threads=params.threads,
        k=params.k,
        rare=params.rare,
    )


if __name__ == "__main__":
    main()
