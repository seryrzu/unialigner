import argparse
from collections import namedtuple, defaultdict, Counter
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as spstats

from utils.bio import parse_cigar
from utils.os_utils import expandpath, smart_makedirs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-l", "--hor-len", type=int)
    params = parser.parse_args()


    for input_file in params.input.split(','):
        cx, cy = 0, 0
        x, y = [cx], [cy]
        stats = defaultdict(list)
        with open(input_file) as f:
            parsed_cigar, cnt = parse_cigar(f.readline())
            for length, mode in parsed_cigar:
                stats[mode].append(length)
        # mode_name = {"M": "Match", "X": "Mismatch", "I": "Insertion", "D": "Deletion"}
        mode_name = {"I": "Insertion", "D": "Deletion"}
        for m, mode in mode_name.items():
            plt.hist(stats[m], bins=500, alpha = 0.8)
        # plt.title(mode)
        plt.xlabel("Length")
        plt.ylabel("Count")
        plt.xlim(0, 20000)
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
            print(spstats.poisson.pmf([1, 2, 3, 4], mu=average - 1, loc=1))



main()