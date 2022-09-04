import argparse
import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import sys
from collections import defaultdict
from itertools import groupby

from standard_logger import get_logger
from utils.git import get_git_revision_short_hash
from utils.os_utils import expandpath, smart_makedirs

SCRIPT_FN = os.path.realpath(__file__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--shortest-matches", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--asm1-name", default="")
    parser.add_argument("--asm2-name", default="")
    parser.add_argument("--min-length", default=20, type=int)
    parser.add_argument("--edlib", default="")
    params = parser.parse_args()

    params.outdir = expandpath(params.outdir)
    smart_makedirs(params.outdir)
    return params


def get_intervals(params, df):
    Xs, Ys = defaultdict(list), defaultdict(list)
    for i, row in df.iterrows():
        if row.Length < params.min_length:
            continue
        freqs = (row.FstFreq, row.SndFreq)
        for Xst in str(row.FstStarts).split(','):
            Xst = int(Xst)
            for Yst in str(row.SndStarts).split(','):
                Yst = int(Yst)
                Xs[freqs].append(Xst)
                Xs[freqs].append(Xst + row.Length)
                Ys[freqs].append(Yst)
                Ys[freqs].append(Yst + row.Length)
                Xs[freqs].append(np.nan)
                Ys[freqs].append(np.nan)
    return Xs, Ys


def parse_cigar(cigar_fn):
    X, Y = [0], [0]
    with open(cigar_fn) as f:
        cigar = f.readline().strip()
    cig_iter = groupby(cigar, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op == 'M' or op == 'X':
            X.append(X[-1] + length)
            Y.append(Y[-1] + length)
        elif op == 'D':
            X.append(X[-1] + length)
            Y.append(Y[-1])
        else:
            assert op == 'I'
            X.append(X[-1])
            Y.append(Y[-1] + length)
    return X, Y


def get_traces(al_X, al_Y, edlib_X, edlib_Y, Xs, Ys):
    traces = [go.Scattergl(x=al_X, y=al_Y, mode='lines', name='Alignment', legendgroup='Alignment',
                           marker=dict(size=2),
                           line=dict(color='black', width=4),
                           showlegend=True)]
    if edlib_Y is not None:
        traces.append(go.Scattergl(x=edlib_X, y=edlib_Y, mode='lines', name='Edlib Alignment', legendgroup='Edlib Alignment',
                                   marker=dict(size=2),
                                   line=dict(color='red', width=4),
                                   showlegend=True))
    freqs = list(Xs.keys())
    for i, freq in enumerate(freqs[::-1]):
        X = Xs[freq]
        Y = Ys[freq]
        traces.append(go.Scattergl(x=X,
                                   y=Y,
                                   mode='markers+lines',
                                   name=str(freq),
                                   legendgroup=str(freq),
                                   marker=dict(size=2),
                                   line=dict(width=1),
                                   showlegend=True))
    return traces


def main():
    params = parse_args()

    logfn = os.path.join(params.outdir, 'minseq.log')
    global logger
    logger = get_logger(logfn,
                        logger_name='minseq')
    logger.info(f'Constructing dotplot with {SCRIPT_FN} started')
    logger.info('cmd: {}'.format(sys.argv))
    logger.info('git hash: {}'.format(get_git_revision_short_hash()))

    df = pd.read_csv(os.path.join(params.shortest_matches, 'shortest_matches.tsv'), header=0, sep='\t')

    Xs, Ys = get_intervals(params, df)
    al_X, al_Y = parse_cigar(os.path.join(params.shortest_matches, 'cigar.txt'))
    has_edlib = params.edlib != ''
    edlib_X, edlib_Y = None, None
    if has_edlib:
        edlib_X, edlib_Y = parse_cigar(os.path.join(params.shortest_matches, params.edlib))

    # n_colors = 1 + len(Xs)
    # colors = px.colors.sample_colorscale("plotly3", [n / (n_colors - 1) for n in range(n_colors)])

    traces = get_traces(al_X, al_Y, edlib_X, edlib_Y, Xs, Ys)

    layout = go.Layout(title='Dotplot', xaxis_title=params.asm1_name, yaxis_title=params.asm2_name)
    fig = go.Figure(traces, layout)
    fig.update_yaxes(scaleanchor = "x", scaleratio = 1)
    fig.write_html(os.path.join(params.outdir, 'dotplot.html'))


if __name__ == "__main__":
    main()
