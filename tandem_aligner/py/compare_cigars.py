import argparse
from collections import namedtuple
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from utils.bio import parse_cigar
from utils.os_utils import expandpath, smart_makedirs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i1", required=True)
    parser.add_argument("-i2", required=True)
    params = parser.parse_args()

    all_coords = []
    for input_file in [params.i1, params.i2]:
        all_coords.append(set())
        matched_coords = all_coords[-1]
        cx, cy = 0, 0
        with open(input_file) as f:
            parsed_cigar, cnt = parse_cigar(f.readline())
            for length, mode in parsed_cigar:
                if mode in "M=X":
                    if mode in "M=":
                        for i in range(length):
                            matched_coords.add((cx + i, cy + i))
                    cx += length
                    cy += length
                elif mode == "I":
                    cy += length
                elif mode == "D":
                    cx += length

    print(len(all_coords[0]), len(all_coords[1]))
    print(len(all_coords[0].intersection(all_coords[1])))
    print(
        len(all_coords[0].intersection(all_coords[1]))
        / min(len(all_coords[0]), len(all_coords[1]))
    )


main()
