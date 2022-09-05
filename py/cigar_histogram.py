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
    for input_file in params.input.split(','):
        stats = defaultdict(list)
        with open(input_file) as f:
            parsed_cigar, cnt = parse_cigar(f.readline())
            for length, mode in parsed_cigar:
                stats[mode].append(length)
        # mode_name = {"M": "Match", "X": "Mismatch", "I": "Insertion", "D": "Deletion"}
        mode_name = {"I": "Insertion", "D": "Deletion"}
        for m, mode in mode_name.items():
            plt.hist(stats[m], bins = np.arange(0, 20000 + 200, 200), alpha = 0.8)
        # plt.title(mode)
        plt.xlabel("Length")
        plt.ylabel("Count")
        # plt.xlim(0, 20000)
        plt.legend(list(mode_name.values()))
        plt.savefig(os.path.join(params.output, input_file + "_hist.pdf"), format="pdf")
        plt.close()
        if params.hor_len:
            hor_len = params.hor_len
            hor_indel_lens = []
            for m, mode in mode_name.items():
                for k in stats[m]:
                    if abs(k - round(k/hor_len) * hor_len) < 0.05:
                        hor_indel_lens.append(int(round(k/hor_len)))

            counter = Counter(hor_indel_lens)
            print("Counter of hor_indel_lens:\n")
            for i in range(1, max(counter) + 1):
                print(i, counter[i])
            print("\nNormalized counter (frequences):")
            for i in range(1, max(counter) + 1):
                print(i, int(round(counter[i] / len(hor_indel_lens) * 100)))

            average = sum(hor_indel_lens) / len(hor_indel_lens)
            print("Average length (estimated parameter for Poisson): ", average)
            print(spstats.poisson.pmf([1, 2, 3, 4, 5, 6], mu=average - 1, loc=1))


def est_mism_shortindel_rate(params):
    for input_file in params.input.split(','):
        with open(input_file) as f:
            parsed_cigar, cnt = parse_cigar(f.readline())

        norm_alignment_len = 0
        len_mism = 0
        len_short_indels= 0
        for length, mode in parsed_cigar:
            if mode == "X":
                len_mism += length
            if mode in "ID":
                if length > params.max_short_indel_len:
                    continue
                len_short_indels += length
            norm_alignment_len += length
        print(norm_alignment_len,
              len_mism,
              len_short_indels,
              len_mism / norm_alignment_len,
              len_short_indels / norm_alignment_len)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-l", "--hor-len", type=int)
    parser.add_argument("-m", "--max-short-indel-len", type=int, default=5)
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
        # plt.savefig(os.path.join(params.output, input_file + "_mism_hist.pdf"), format="pdf")
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
        # plt.savefig(os.path.join(params.output, input_file + "_smallindel_hist.pdf"), format="pdf")
        # plt.close()

main()