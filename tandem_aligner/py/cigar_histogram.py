import argparse
from collections import namedtuple, defaultdict, Counter
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as spstats

from utils.bio import parse_cigar
from utils.os_utils import expandpath, smart_makedirs


def plot_hist(params):
    for input_file in params.input.split(","):
        stats = defaultdict(list)
        with open(input_file) as f:
            parsed_cigar, cnt = parse_cigar(f.readline())
            i = 0
            for length, mode in parsed_cigar:
                stats[mode].append((length, i))
                if mode != "I":
                    i += length
        # mode_name = {"M": "Match", "X": "Mismatch", "I": "Insertion", "D": "Deletion"}
        mode_name = {"I": "Insertion", "D": "Deletion"}
        max_len = max([len for len, pos in stats[mode] for mode in mode_name])
        for m, mode in mode_name.items():
            print(
                f"{len(stats[m])} {mode} blocks account for {sum([len for len,pos in stats[m]])} {mode.lower()}s"
            )
            plt.hist(
                [x[0] for x in stats[m]],
                bins=np.arange(0, max_len + params.bin_size, params.bin_size),
                alpha=0.5,
            )
        # plt.title(mode)
        plt.xlabel("Length (bp)")
        plt.ylabel("Count")
        # plt.xlim(0, 20000)
        plt.legend(list(mode_name.values()))
        if params.ylog:
            plt.yscale("log")
        if params.xlog:
            plt.xscale("log")
        plt.savefig(
            os.path.join(params.outdir, input_file.split("/")[-1][:-3] + "_hist.pdf"),
            format="pdf",
        )
        plt.close()
        if params.hor_len:
            hor_len = params.hor_len
            hor_indel_lens = []
            norm2real = defaultdict(list)
            hor_indel2pos = defaultdict(list)
            for m, mode in mode_name.items():
                for k, pos in stats[m]:
                    if (
                        round(k / hor_len)
                        and abs(k - round(k / hor_len) * hor_len) / hor_len < 0.05
                    ):
                        hor_mult = int(round(k / hor_len))
                        hor_indel_lens.append(hor_mult)
                        norm2real[hor_mult].append(k)
                        hor_indel2pos[hor_mult].append(pos)

            counter = Counter(hor_indel_lens)
            print("Counter of hor_indel_lens:\n")
            for i in range(1, max(counter) + 1):
                print(i, counter[i], sorted(norm2real[i]))
            print("\nNormalized counter (frequences):")
            for i in range(1, max(counter) + 1):
                print(i, int(round(counter[i] / len(hor_indel_lens) * 100)))

            average = sum(hor_indel_lens) / len(hor_indel_lens)
            print("Average length (estimated parameter for Poisson): ", average)
            print(spstats.poisson.pmf([1, 2, 3, 4, 5, 6], mu=average - 1, loc=1))

            for height, positions in sorted(hor_indel2pos.items()):
                plt.scatter(
                    np.array(positions) / 1e6, [height] * len(positions), marker="."
                )
            plt.legend(sorted(hor_indel2pos))
            plt.ylabel("HOR-indel multiplicity")
            plt.xlabel("cenX-1 (Mb)")
            plt.savefig(
                os.path.join(
                    params.outdir, input_file.split("/")[-1][:-3] + "_hor_indels.pdf"
                ),
                format="pdf",
            )


def est_mism_shortindel_rate(params):
    for input_file in params.input.split(","):
        with open(input_file) as f:
            parsed_cigar, cnt = parse_cigar(f.readline())

        norm_alignment_len = 0
        (
            n_mism,
            n_short_insertions,
            n_short_deletions,
            n_long_insertions,
            n_long_deletions,
        ) = (0, 0, 0, 0, 0)
        (
            len_mism,
            len_short_insertions,
            len_short_deletions,
            len_long_insertions,
            len_long_deletions,
        ) = (0, 0, 0, 0, 0)
        max_insertion, max_deletion = 0, 0
        for length, mode in parsed_cigar:
            if mode == "X":
                n_mism += 1
                len_mism += length
            if mode == "I":
                if length > params.max_short_indel_len:
                    n_long_insertions += 1
                    len_long_insertions += length
                    max_insertion = max(max_insertion, length)
                    continue
                else:
                    n_short_insertions += 1
                    len_short_insertions += length
            if mode == "D":
                if length > params.max_short_indel_len:
                    n_long_deletions += 1
                    len_long_deletions += length
                    max_deletion = max(max_deletion, length)
                    continue
                else:
                    n_short_deletions += 1
                    len_short_deletions += length
            norm_alignment_len += length

        print("norm_alignment_len = ", norm_alignment_len)
        print("n_mism = ", n_mism, "len_mism = ", len_mism)
        print(
            "n_short_insertions = ",
            n_short_insertions,
            "len_short_indels = ",
            len_short_insertions,
        )
        print(
            "n_long_insertions = ",
            n_long_insertions,
            "len_long_indels = ",
            len_long_insertions,
        )
        print(
            "n_short_deletions = ",
            n_short_deletions,
            "len_short_indels = ",
            len_short_deletions,
        )
        print(
            "n_long_deletions = ",
            n_long_deletions,
            "len_long_indels = ",
            len_long_deletions,
        )
        print("len_mism / norm_alignment_len", len_mism / norm_alignment_len)
        print(
            "len_short_insertions / norm_alignment_len",
            len_short_insertions / norm_alignment_len,
        )
        print("len_long_insertions", len_long_insertions)
        print(
            "len_short_deletions / norm_alignment_len",
            len_short_deletions / norm_alignment_len,
        )
        print("len_long_deletions", len_long_deletions)
        print("max_insertion", max_insertion)
        print("max_deletion", max_deletion)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--outdir")
    parser.add_argument("-l", "--hor-len", type=int)
    parser.add_argument("-m", "--max-short-indel-len", type=int, default=5)
    parser.add_argument("--xlog", action="store_true")
    parser.add_argument("--ylog", action="store_true")
    parser.add_argument("--bin_size", type=int, default=200)
    params = parser.parse_args()

    plot_hist(params)
    est_mism_shortindel_rate(params)

    # pos_mism = []
    # pos = 0
    # for length, mode in parsed_cigar:
    #     if mode == "X":
    #         pos_mism.append(pos)
    #         pos += length
    #     if mode == "M" or mode == "D":
    #         pos += length

    # plt.hist(pos_mism, bins=100)
    # plt.savefig(os.path.join(params.outdir, input_file + "_mism_hist.pdf"), format="pdf")
    # plt.close()

    # pos_short_indel = []
    # pos = 0
    # for length, mode in parsed_cigar:
    #     if mode in "DI":
    #         if length <= params.max_small_indel_len:
    #             pos_short_indel.append(pos)
    #     if mode in "MDX":
    #         pos += length

    # plt.hist(pos_short_indel, bins=100)
    # plt.savefig(os.path.join(params.outdir, input_file + "_smallindel_hist.pdf"), format="pdf")
    # plt.close()


main()
