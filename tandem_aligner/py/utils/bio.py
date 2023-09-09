# (c) 2020 by Authors
# This file is a part of the SD program.
# see LICENSE file

from collections import OrderedDict
from itertools import groupby

from Bio import SeqIO
import re


def read_bio_seq(filename):
    seqs = read_bio_seqs(filename)
    return str(list(seqs.values())[0])


def read_bio_seqs(filename):
    form = filename.split(".")[-1]
    if form == "fa" or form == "fna":
        form = "fasta"
    elif form == "fq":
        form = "fastq"
    raw_seqs = SeqIO.parse(filename, format=form)
    seqs = OrderedDict()
    for seq in raw_seqs:
        seqs[seq.id] = str(seq.seq)
    return seqs


trans = str.maketrans("ATGCatgc-", "TACGtacg-")


def RC(s):
    return s.translate(trans)[::-1]


def write_bio_seqs(filename, seqs, width=60):
    with open(filename, "w") as f:
        for seq_id, seq in seqs.items():
            print(">{}".format(seq_id), file=f)
            if width is None:
                print(seq, file=f)
                continue
            start = 0
            seq_len = len(seq)
            while seq_len - start >= width:
                print(seq[start : start + width], file=f)
                start += width
            print(seq[start:], file=f)


def parse_cigar(cigar, s1=None, s2=None):
    parsed_cigar = []
    st = 0
    cnt = dict.fromkeys(list("=MXID"), 0)
    for mo in re.finditer(r"=|M|X|I|D", cigar):
        group = mo.group()
        pos = mo.start()
        region_len = int(cigar[st:pos])
        parsed_cigar.append((region_len, group))
        cnt[group] += region_len
        st = pos + 1
    if s1 is None or s2 is None:
        return parsed_cigar, cnt

    a1, a2 = [], []
    i1, i2 = 0, 0
    for region_len, group in parsed_cigar:
        if group in "=XM":
            new_s1 = s1[i1 : i1 + region_len]
            new_s2 = s2[i2 : i2 + region_len]
            # if group == '=':
            #     assert new_s1 == new_s2
            a1 += new_s1
            a2 += new_s2
            i1 += region_len
            i2 += region_len
        elif group == "D":
            a1 += "-" * region_len
            a2 += s2[i2 : i2 + region_len]
            i2 += region_len
        elif group == "I":
            a2 += "-" * region_len
            a1 += s1[i1 : i1 + region_len]
            i1 += region_len

    # a1 = ''.join(a1)
    # a2 = ''.join(a2)
    return parsed_cigar, cnt, a1, a2


assert parse_cigar("89=1X6=3X76=") == (
    [(89, "="), (1, "X"), (6, "="), (3, "X"), (76, "=")],
    {"=": 171, "X": 4, "I": 0, "D": 0, "M": 0},
)


def calc_identity(a, b, mode="NW", k=None):
    import edlib

    alignment = edlib.align(a, b, task="path", mode=mode, k=k)
    if alignment["editDistance"] == -1:
        return 0, alignment
    cigar, cigar_stats = parse_cigar(alignment["cigar"])
    alignment_len = sum(cigar_stats.values())
    identity = 1 - alignment["editDistance"] / alignment_len
    assert 0 <= identity <= 1
    return identity, alignment


def compress_homopolymer(seq, return_list=False):
    compressed_seq = [x[0] for x in groupby(list(seq))]
    if return_list:
        return compressed_seq
    return "".join(compressed_seq)
