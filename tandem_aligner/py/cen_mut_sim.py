import argparse
import numpy as np
from collections import namedtuple
from utils.os_utils import smart_makedirs
from utils.bio import read_bio_seq, write_bio_seqs
from itertools import groupby
import os


class HORBlock:
    def __init__(self, line, seq):
        self.ref_name, self.monomers, self.st, en, self.ident = line.strip().split()
        self.st = int(self.st)
        self.ident = float(self.ident)


def def_lam():
    return 1.6594594594594594


def read_hors(hor_fn, seq_fn):
    seq = read_bio_seq(seq_fn)
    hors = []
    with open(hor_fn) as f:
        for line in f:
            hors.append(HORBlock(line, seq))
    return hors


def mutate(seq, mism=0.001, delet=0.0001, ins=0.0001):
    mutseq = []
    uncompr_cigar = []
    choice = np.random.choice(
        list("mxdi"), size=len(seq), p=[1 - mism - delet - ins, mism, delet, ins]
    )
    for char, ch in zip(seq, choice):
        if ch == "m":
            uncompr_cigar.append("M")
            mutseq.append(char)
        elif ch == "x":
            uncompr_cigar.append("X")
            bases = list("ACGT")
            bases.remove(char)
            mutseq.append(np.random.choice(bases, size=1)[0])
        elif ch == "d":
            uncompr_cigar.append("D")
        else:
            bases = list("ACGT")
            mutseq.append(np.random.choice(bases, size=1)[0])
            uncompr_cigar.append("I")
            mutseq.append(char)
            uncompr_cigar.append("M")
    assert len(mutseq) == sum(x in "MXI" for x in uncompr_cigar)
    assert len(seq) == sum(x in "MXD" for x in uncompr_cigar)
    return "".join(mutseq), uncompr_cigar


def introduce_cnv(hor_string, hors, seq_fn, n_cnv, lam, dup_ratio):
    seq = read_bio_seq(seq_fn)
    hor_string_length = len(hor_string)
    lens = 1 + np.random.poisson(lam=lam - 1, size=n_cnv)
    for length in lens:
        st = np.random.randint(0, len(hor_string) - length - 1, size=1)[0]
        is_dup = np.random.random(1)[0] < dup_ratio
        if is_dup:
            hor_string = hor_string[: st + length] + hor_string[st:]
        else:
            hor_string = hor_string[:st] + hor_string[st + length :]

    uncomp_cigar = ["M"] * hors[0].st
    cur_max = -1
    MockHOR = namedtuple("MockHOR", ["st"])
    hors.append(MockHOR(len(seq)))
    mod_seq = []
    for i in hor_string:
        print(i)
        subseq = seq[hors[i].st : hors[i + 1].st]
        mut_subseq, mut_subcigar = mutate(subseq)
        mod_seq.append(mut_subseq)
        if cur_max + 1 == i:
            uncomp_cigar += mut_subcigar
            cur_max = i
        elif cur_max >= i:
            uncomp_cigar += len(mut_subseq) * "I"
        else:
            for j in range(cur_max + 1, i):
                uncomp_cigar += (hors[j + 1].st - hors[j].st) * "D"
            uncomp_cigar += mut_subcigar
            cur_max = i

    assert len(seq) == sum(x in "MXD" for x in uncomp_cigar)

    mod_seq = "".join(mod_seq)
    assert len(mod_seq) == sum(x in "MXI" for x in uncomp_cigar)

    comp_cigar = []
    for k, g in groupby(uncomp_cigar):
        comp_cigar.append((k, len(list(g))))

    cigar = "".join(str(v) + str(k) for k, v in comp_cigar)
    return mod_seq, cigar


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--seq", required=True)
    parser.add_argument("-d", "--hors", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-N", "--n-cnv", type=int, required=True)
    parser.add_argument("--lam", type=float, default=def_lam())
    parser.add_argument("--dup", type=float, default=0.5)
    params = parser.parse_args()

    smart_makedirs(params.outdir)
    np.random.seed(params.seed)
    hors = read_hors(params.hors, params.seq)
    hors[0].st = 0
    hor_string = list(range(len(hors)))

    mod_seq, cigar = introduce_cnv(
        hor_string, hors, params.seq, params.n_cnv, params.lam, params.dup
    )

    with open(os.path.join(params.outdir, "true_cigar.txt"), "w") as f:
        print(cigar, file=f)

    write_bio_seqs(os.path.join(params.outdir, "mod.fasta"), {"mod": mod_seq})


if __name__ == "__main__":
    main()
