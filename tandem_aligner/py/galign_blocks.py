import argparse
from collections import namedtuple
import os
import subprocess
import sys

from standard_logger import get_logger
from utils.git import get_git_revision_short_hash
from utils.bio import parse_cigar
from utils.os_utils import smart_makedirs


ComprBlock = namedtuple("ComprBlock", ["len1", "len2", "status"])
AlignedBlock = namedtuple("AlignedBlock", ["st1", "en1", "st2", "en2"])


def get_compr_cigar(parsed_cigar):
    compr_cigar = []
    i = 0
    while i < len(parsed_cigar):
        _, block_status = parsed_cigar[i]
        len1, len2 = 0, 0
        while i < len(parsed_cigar):
            length, status = parsed_cigar[i]
            if block_status == status == "=":
                len1 += length
                len2 += length
                i += 1
            elif block_status != "=" and status != "=":
                if status == "X":
                    len1 += length
                    len2 += length
                elif status == "D":
                    len2 += length
                else:
                    assert status == "I"
                    len1 += length
                i += 1
            else:
                break
        block_status = "=" if block_status == "=" else "X"
        compr_cigar.append(ComprBlock(len1, len2, block_status))
    return compr_cigar


def cigar2blocks(cigar, tol_gap, logger, outfn=None, min_large_unmatched_block=10):
    parsed_cigar, _ = parse_cigar(cigar)

    compr_cigar = get_compr_cigar(parsed_cigar)
    logger.info("Large (size > {min_large_unmatched_block}) unmatched blocks:")
    for block in compr_cigar:
        if (
            block.status != "="
            and max(block.len1, block.len2) > min_large_unmatched_block
        ):
            logger.info(block)
    logger.info("End large unmatched blocks:")

    i = 0
    st1, st2 = 0, 0
    if compr_cigar[0].status != "=":
        st1 += compr_cigar[0].len1
        st2 += compr_cigar[0].len2
        i += 1

    if len(compr_cigar[i:]) % 2 == 1:
        compr_cigar.append(ComprBlock(0, 0, "X"))

    assert compr_cigar[i].status == "="

    blocks = []
    en1, en2 = st1, st2
    for bl1, bl2 in zip(compr_cigar[i::2], compr_cigar[i + 1 :: 2]):
        last_add = False
        en1 += bl1.len1
        en2 += bl1.len2
        if max(bl2.len1, bl2.len2) > tol_gap:
            blocks.append(AlignedBlock(st1, en1, st2, en2))
            en1 += bl2.len1
            en2 += bl2.len2
            st1, st2 = en1, en2
            last_add = True
        else:
            en1 += bl2.len1
            en2 += bl2.len2
    if not last_add:
        blocks.append(AlignedBlock(st1, en1, st2, en2))

    if outfn is not None:
        with open(outfn, "w") as f:
            print(f"st1\ten1\tst2\ten2", file=f)
            for block in blocks:
                print(block.st1, block.en1, block.st2, block.en2, file=f, sep="\t")
    return blocks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-1", "--s1", required=True)
    parser.add_argument("-2", "--s2", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--tol-gap", required=True, type=int)
    params = parser.parse_args()

    smart_makedirs(params.outdir)

    logfn = os.path.join(params.outdir, "galign_blocks.log")
    logger = get_logger(logfn, logger_name="galign_blocks")
    SCRIPT_FN = os.path.realpath(__file__)
    logger.info(f"Aligning with {SCRIPT_FN} started")
    logger.info("cmd: {}".format(sys.argv))
    logger.info("git hash: {}".format(get_git_revision_short_hash()))

    edlib_cmd = f"edlib-aligner -p -f CIG_EXT {params.s1} {params.s2}".split()
    edlib_proc = subprocess.Popen(edlib_cmd, stdout=subprocess.PIPE)
    edlib_output = edlib_proc.stdout.read().decode("utf-8").split("\n")
    logger.info("Edlib results")
    for line in edlib_output:
        logger.info(line)

    cigar = edlib_output[-4]
    with open(os.path.join(params.outdir, "cigar.txt"), "w") as f:
        print(cigar, file=f)

    cigar2blocks(
        cigar,
        params.tol_gap,
        logger=logger,
        outfn=os.path.join(params.outdir, "blocks.tsv"),
    )


if __name__ == "__main__":
    main()
