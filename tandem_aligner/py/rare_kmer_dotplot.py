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
    parser.add_argument("--alignment", default="")
    parser.add_argument("--axis-mult", type=float, default=1.0)
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
    if len(params.alignment):
        params.alignment = expandpath(params.alignment)

    smart_makedirs(params.outdir)
    return params


def get_kmer2coords(s, kmers):
    k = len(next(iter(kmers)))
    kmer2coords = defaultdict(list)
    for i in range(len(s) - k + 1):
        kmer = s[i : i + k]
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


def get_alignment(filepath, k, asm1, asm2):
    if len(filepath) == 0:
        return None, None
    raw_Xs, raw_Ys = [0], [0]
    with open(filepath) as f:
        for line in f:
            x, y = line.strip().split()
            if x != "-" and y != "-":
                x, y = map(int, (x, y))
                raw_Xs.append(x)
                raw_Ys.append(y)

    raw_Xs.append(len(asm1))
    raw_Ys.append(len(asm2))
    Xs, Ys = [0], [0]
    for x, y in zip(raw_Xs, raw_Ys):
        K = k
        while K and (Xs[-1] + K > x or Ys[-1] + K > y):
            K -= 1

        while (
            Xs[-1] + K < x and Ys[-1] + K < y and asm1[Xs[-1] + K] == asm2[Ys[-1] + K]
        ):
            K += 1

        while x > Xs[-1] + K and y > Ys[-1] + K and asm1[x - 1] == asm2[y - 1]:
            x -= 1
            y -= 1
        assert Xs[-1] + K <= x
        assert Ys[-1] + K <= y

        Xs.append(Xs[-1] + K)
        Ys.append(Ys[-1] + K)
        Xs.append(x)
        Ys.append(Ys[-1])
        Xs.append(x)
        Ys.append(y)

    filt_X, filt_Y = [Xs[0]], [Ys[0]]
    for x, y, nx, ny in zip(Xs[1:-1], Ys[1:-1], Xs[2:], Ys[2:]):
        if x - y != nx - ny:
            filt_X.append(x)
            filt_Y.append(y)
            filt_X.append(nx)
            filt_Y.append(ny)
    filt_X.append(Xs[-1])
    filt_Y.append(Ys[-1])

    return filt_X, filt_Y


def main():
    params = parse_args()

    logfn = os.path.join(params.outdir, "dotplot.log")
    global logger
    logger = get_logger(logfn, logger_name="dotplot")
    logger.info(f"Constructing dotplot with {SCRIPT_FN} started")
    logger.info("cmd: {}".format(sys.argv))
    logger.info("git hash: {}".format(get_git_revision_short_hash()))

    asm1 = read_bio_seqs(params.asm1)
    asm2 = read_bio_seqs(params.asm2)

    asm1_name, asm1 = next(iter(asm1.items()))
    asm2_name, asm2 = next(iter(asm2.items()))
    asm1_name, asm2_name = asm1_name.split("|")[0], asm2_name.split("|")[0]

    if params.asm1_name != "":
        asm1_name = params.asm1_name
    if params.asm2_name != "":
        asm2_name = params.asm2_name

    X_al, Y_al = get_alignment(params.alignment, params.k, asm1, asm2)

    rare_kmers1 = get_rare_kmers(asm1, params.k, params.rare)
    rare_kmers2 = get_rare_kmers(asm2, params.k, params.rare)
    rare_kmers = rare_kmers1 & rare_kmers2

    kmer2coords1 = get_kmer2coords(asm1, rare_kmers)
    kmer2coords2 = get_kmer2coords(asm2, rare_kmers)

    pairs = get_pairs(kmer2coords1, kmer2coords2)
    plt.figure(figsize=(5, 5), dpi=200)
    plt.plot(
        [pair[0] * params.axis_mult for pair in pairs],
        [pair[1] * params.axis_mult for pair in pairs],
        marker=".",
        ms=0.5,
        alpha=1,
        color="black",
    )
    if params.alignment != "":
        plt.plot(
            [x * params.axis_mult for x in X_al],
            [y * params.axis_mult for y in Y_al],
            color="red",
        )
    plt.xlabel(asm1_name)
    plt.ylabel(asm2_name)
    plt.xlim(0)
    plt.ylim(0)
    plt.grid(b=True, which="major", axis="both")
    plt.xticks(rotation=30)
    plt.title(f"k={params.k}, MaxCount={params.rare}")
    outfn = os.path.join(params.outdir, "dotplot.png")
    plt.savefig(outfn, format="png", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
