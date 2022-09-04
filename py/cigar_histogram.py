import argparse
from collections import namedtuple, defaultdict
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from utils.bio import parse_cigar
from utils.os_utils import expandpath, smart_makedirs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    params = parser.parse_args()


    for input_file in params.input.split(','):
        cx, cy = 0, 0
        x, y = [cx], [cy]
        stats = defaultdict(list)
        with open(input_file) as f:
            parsed_cigar, cnt = parse_cigar(f.readline())
            for length, mode in parsed_cigar:
                stats[mode].append(length)
        print(stats)
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


main()