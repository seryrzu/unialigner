import argparse
import numpy as np
from collections import namedtuple
from utils.os_utils import smart_makedirs
from utils.bio import read_bio_seq, write_bio_seqs
from itertools import groupby
import os

from cen_mut_sim import mutate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--seq", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-m", "--mut", type=float, default=0.02)
    parser.add_argument("-d", "--del-len", type=int, default=1000)
    params = parser.parse_args()

    smart_makedirs(params.outdir)
    np.random.seed(params.seed)
    seq = read_bio_seq(params.seq)

    del_pos = np.random.randint(0, len(seq) - params.del_len, 1)[0]
    prefix, suffix = seq[:del_pos], seq[del_pos + params.del_len :]
    mut = params.mut
    mut_prefix, uncompr_cigar_prefix = mutate(
        prefix, mism=mut / 2, delet=mut / 4, ins=mut / 4
    )
    mut_suffix, uncompr_cigar_suffix = mutate(
        suffix, mism=mut / 2, delet=mut / 4, ins=mut / 4
    )

    uncompr_cigar = uncompr_cigar_prefix + ["D"] * params.del_len + uncompr_cigar_suffix
    mut_seq = mut_prefix + mut_suffix

    cigar = []
    for k, g in groupby(uncompr_cigar):
        cigar.append((k, len(list(g))))
    cigar = "".join(str(v) + str(k) for k, v in cigar)

    with open(os.path.join(params.outdir, "true_cigar.txt"), "w") as f:
        print(cigar, file=f)

    write_bio_seqs(os.path.join(params.outdir, "mod.fasta"), {"mod": mut_seq})


if __name__ == "__main__":
    main()
