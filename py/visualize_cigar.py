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
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("--asm1-name", default="")
    parser.add_argument("--asm2-name", default="")
    parser.add_argument("--legend", default="")
    params = parser.parse_args()


    for input_file in params.input.split(','):
        cx, cy = 0, 0
        x, y = [cx], [cy]
        with open(input_file) as f:
            parsed_cigar = parse_cigar(f.readline())[0]
            for len, mode in parsed_cigar:
                len = int(len)
                if mode == 'M' or mode == 'S':
                    cx += len
                    cy += len
                elif mode == 'I':
                    cy += len
                elif mode == 'D':
                    cx += len

                x.append(cx)
                y.append(cy)

        plt.plot(x, y)

    if params.legend:
        plt.legend(params.legend.split(','))
    plt.xlabel(params.asm1_name)
    plt.ylabel(params.asm2_name)
    plt.savefig(os.path.join(params.output, 'cigar.pdf'), format='pdf')

main()