from collections import defaultdict
import argparse
import os
import sys
import pandas as pd
import numpy as np

import plotly.graph_objs as go
import plotly.express as px

from standard_logger import get_logger
from utils.os_utils import expandpath, smart_makedirs
from utils.git import get_git_revision_short_hash


SCRIPT_FN = os.path.realpath(__file__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--shortest-matches", required=True)
    parser.add_argument("--asm1-name", default="")
    parser.add_argument("--asm2-name", default="")
    parser.add_argument("-o", "--outdir", required=True)
    params = parser.parse_args()

    params.outdir = expandpath(params.outdir)
    smart_makedirs(params.outdir)
    return params


def main():
    params = parse_args()

    logfn = os.path.join(params.outdir, 'dotplot.log')
    global logger
    logger = get_logger(logfn,
                        logger_name='dotplot')
    logger.info(f'Constructing dotplot with {SCRIPT_FN} started')
    logger.info('cmd: {}'.format(sys.argv))
    logger.info('git hash: {}'.format(get_git_revision_short_hash()))

    asm1_name, asm2_name = "", ""
    if params.asm1_name != "":
        asm1_name = params.asm1_name
    if params.asm2_name != "":
        asm2_name = params.asm2_name

    df = pd.read_csv(params.shortest_matches, header=0, sep='\t')

    Xs, Ys = defaultdict(list), defaultdict(list)
    for i, row in df.iterrows():
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

    layout = go.Layout(
        title='Dotplot',
        xaxis_title=asm1_name,
        yaxis_title=asm2_name
    )

    n_colors = len(Xs)
    colors = px.colors.sample_colorscale("turbo", [n/(n_colors - 1) * 0.8 for n in range(n_colors)])

    traces = []
    freqs = list(Xs.keys())
    for i, freq in enumerate(freqs[::-1]):
        X = Xs[freq]
        Y = Ys[freq]
        traces.append(go.Scattergl(x=X,
                                   y=Y,
                                   mode='markers+lines',
                                   name=str(freq),
                                   legendgroup=str(freq),
                                   line=dict(color=colors[i]),
                                   showlegend=True))
    fig = go.Figure(traces, layout)

    # Plot the chart
    # fig = go.Figure([data], layout)

    # pyo.plot(fig, filename=os.path.join(params.outdir, 'dotplot.html'), auto_open=False)
    fig.write_html(os.path.join(params.outdir, 'dotplot.html'))
    # fig.write_image(os.path.join(params.outdir, 'dotplot.pdf'))


if __name__ == "__main__":
    main()