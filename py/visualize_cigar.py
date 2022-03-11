import argparse
from collections import namedtuple
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from utils.os_utils import expandpath, smart_makedirs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    params = parser.parse_args()

    cx, cy = 0, 0
    x, y = [cx], [cy]
    with open(params.input) as f:
        for line in f:
            mode, len = line.split()
            len = int(len)
            if mode == 'M' or mode == 'S':
                cx += len
                cy += len
            elif mode == 'I':
                cx += len
            elif mode == 'D':
                cy += len

            x.append(cx)
            y.append(cy)

    plt.plot(x, y, 'bo-')
    plt.savefig(os.path.join(params.output, 'cigar.pdf'), format='pdf')

main()