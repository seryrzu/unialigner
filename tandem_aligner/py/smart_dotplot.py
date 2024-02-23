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
from get_bridges import get_bridges
from get_genes_orthologues import get_first_genes, get_second_genes
from get_conserved_regions import get_conserved_traces

SCRIPT_FN = os.path.realpath(__file__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--bridges", action="store_true")
    parser.add_argument("--asm1-name", default="")
    parser.add_argument("--asm2-name", default="")
    parser.add_argument("--min-lengths", default="1,20,50,80")
    parser.add_argument("--edlib", default="")
    parser.add_argument("--ground-truth", default="")
    parser.add_argument("--first-genes", default="")
    parser.add_argument("--second-genes", default="")
    params = parser.parse_args()

    params.outdir = expandpath(params.outdir)
    smart_makedirs(params.outdir)
    return params


def get_intervals(params, df):
    Xs, Ys = defaultdict(list), defaultdict(list)
    visible = defaultdict(lambda: defaultdict(list))

    min_lengths = sorted(
        [int(min_length) for min_length in params.min_lengths.split(",")]
    )
    min_min_length = min_lengths[0]  # minimum min_length

    for i, row in df.iterrows():
        if row.Length < min_min_length:
            continue

        # maximum min_length under interval length
        max_min_length = min_min_length
        for min_length in min_lengths:
            if row.Length > min_length and max_min_length < min_length:
                max_min_length = min_length
        # group = (freq1,freq2,maximum min_length)
        group = (row.FstFreq, row.SndFreq, max_min_length)

        for Xst in str(row.FstStarts).split(","):
            Xst = int(Xst)
            for Yst in str(row.SndStarts).split(","):
                Yst = int(Yst)
                Xs[group].append(Xst)
                Xs[group].append(Xst + row.Length)
                Ys[group].append(Yst)
                Ys[group].append(Yst + row.Length)
                Xs[group].append(np.nan)
                Ys[group].append(np.nan)

    return Xs, Ys


def parse_cigar(cigar_fn):
    X, Y = [0], [0]
    with open(cigar_fn) as f:
        cigar = f.readline().strip()
    cig_iter = groupby(cigar, lambda chr: chr.isdigit())
    segment_lengths = {"M": 0, "X": 0, "I": 0, "D": 0}
    for _, length_digits in cig_iter:
        length = int("".join(length_digits))
        op = next(next(cig_iter)[1])
        if op == "M" or op == "X":
            X.append(X[-1] + length)
            Y.append(Y[-1] + length)
        elif op == "D":
            X.append(X[-1] + length)
            Y.append(Y[-1])
        else:
            assert op == "I"
            X.append(X[-1])
            Y.append(Y[-1] + length)
        segment_lengths[op] += length
    print(segment_lengths)
    return X, Y


def get_traces(
    al_X,
    al_Y,
    edlib_X,
    edlib_Y,
    gt_X,
    gt_Y,
    Xs,
    Ys,
    bridges_X,
    bridges_Y,
    first_genes,
    second_genes,
    conserved_traces,
):
    traces = []
    groups = list(Xs.keys())
    groups = sorted(groups, key=lambda x: x[2])
    interval_sizes = sorted({group[2] for group in groups})
    colors = ["yellow", "orange", "red", "brown"]
    color_map = dict()
    for i in range(len(interval_sizes)):
        color_map[interval_sizes[i]] = colors[min(i, len(colors) - 1)]
    for i, group in enumerate(groups):
        X = Xs[group]
        Y = Ys[group]

        group_name = f"{str(group[:2])}-uniqkmer, size:{group[2]}-{'...' if (group[2]==interval_sizes[-1]) or (i==len(groups)-1) else str(groups[i+1][2])}"
        # only first interval adds to legend
        traces.append(
            go.Scattergl(
                x=X,
                y=Y,
                mode="markers+lines",
                name=group_name,
                legendgroup=group_name,
                marker=dict(size=4, color=color_map[group[2]]),
                line=dict(width=3, color=color_map[group[2]]),
                showlegend=True,
            )
        )

    traces.append(
        go.Scattergl(
            x=al_X,
            y=al_Y,
            mode="lines",
            name="TandemAligner",
            legendgroup="TandemAligner",
            marker=dict(size=1),
            line=dict(color="black", width=2),
            showlegend=True,
        )
    )
    if edlib_Y is not None:
        traces.append(
            go.Scattergl(
                x=edlib_X,
                y=edlib_Y,
                mode="lines",
                name="Edlib",
                legendgroup="Edlib",
                marker=dict(size=2),
                line=dict(color="red", width=2),
                showlegend=True,
            )
        )
    if gt_Y is not None:
        traces.append(
            go.Scattergl(
                x=gt_X,
                y=gt_Y,
                mode="lines",
                name="Ground Truth",
                legendgroup="Ground Truth",
                marker=dict(size=3),
                line=dict(color="green", width=2),
                showlegend=True,
            )
        )
    if bridges_X is not None:
        traces.append(
            go.Scattergl(
                x=bridges_X,
                y=bridges_Y,
                mode="markers+lines",
                name="bridges",
                legendgroup="bridges",
                marker=dict(size=6),
                line=dict(width=4, color="pink"),
                showlegend=True,
            )
        )

    if first_genes is not None:
        showlegend = True  # only first iteration will be true
        for (
            first_geneX,
            first_geneY,
            first_orthoX,
            first_orthoY,
            name,
            color,
        ) in first_genes:
            traces.append(
                go.Scattergl(
                    x=first_geneX,
                    y=first_geneY,
                    mode="lines",
                    name="x axis genes",
                    legendgroup="x axis genes",
                    marker=dict(size=1),
                    line=dict(color=color, width=1),
                    showlegend=showlegend,
                    hoverinfo="text",
                    hovertext=name,
                )
            )
            traces.append(
                go.Scattergl(
                    x=first_orthoX,
                    y=first_orthoY,
                    mode="lines",
                    name="projection of x axis genes",
                    legendgroup="projection of x axis genes",
                    marker=dict(size=1),
                    line=dict(color=color, width=1),
                    showlegend=showlegend,
                    hoverinfo="text",
                    hovertext=name,
                    visible="legendonly",
                )
            )
            showlegend = False

    if second_genes is not None:
        showlegend = True  # only first iteration will be true
        for (
            second_geneX,
            second_geneY,
            second_orthoX,
            second_orthoY,
            name,
            color,
        ) in second_genes:
            traces.append(
                go.Scattergl(
                    x=second_geneX,
                    y=second_geneY,
                    mode="lines",
                    name="y axis genes",
                    legendgroup="y axis genes",
                    marker=dict(size=1),
                    line=dict(color=color, width=1),
                    showlegend=showlegend,
                    hoverinfo="text",
                    hovertext=name,
                )
            )
            traces.append(
                go.Scattergl(
                    x=second_orthoX,
                    y=second_orthoY,
                    mode="lines",
                    name="projection of y axis genes",
                    legendgroup="projection of y axis genes",
                    marker=dict(size=1),
                    line=dict(color=color, width=1),
                    showlegend=showlegend,
                    hoverinfo="text",
                    hovertext=name,
                    visible="legendonly",
                )
            )
            showlegend = False

    if conserved_traces is not None:
        traces.append(
            go.Scattergl(
                x=conserved_traces[:, 0],
                y=conserved_traces[:, 1],
                mode="markers+lines",
                name="conserved_regions",
                legendgroup="conserved_regions",
                marker=dict(size=4, color="lime"),
                line=dict(width=3, color="lime"),
                showlegend=True,
            )
        )
    return traces


def main():
    params = parse_args()

    logfn = os.path.join(params.outdir, "tandem_aligner_smart_dotplot.log")
    global logger
    logger = get_logger(logfn, logger_name="tandem_aligner")
    logger.info(f"Constructing dotplot with {SCRIPT_FN} started")
    logger.info("cmd: {}".format(sys.argv))
    logger.info("git hash: {}".format(get_git_revision_short_hash()))

    df = pd.read_csv(
        os.path.join(
            params.indir,
            "shortest_matches.tsv",
        ),
        header=0,
        sep="\t",
    )
    Xs, Ys = get_intervals(params, df)

    al_X, al_Y = parse_cigar(os.path.join(params.indir, "cigar.txt"))
    has_edlib = params.edlib != ""
    edlib_X, edlib_Y = None, None
    if has_edlib:
        edlib_X, edlib_Y = parse_cigar(os.path.join(params.indir, params.edlib))

    gt_X, gt_Y = None, None
    if params.ground_truth:
        gt_X, gt_Y = parse_cigar(os.path.join(params.indir, params.ground_truth))

    bridges_X, bridges_Y = None, None
    if params.bridges:
        bridges_X, bridges_Y = get_bridges(os.path.join(params.indir, "bridges.txt"))

    first_genes = None
    if params.first_genes:
        first_genes = get_first_genes(params.first_genes, al_X, al_Y)

    second_genes = None
    if params.second_genes:
        second_genes = get_second_genes(params.second_genes, al_X, al_Y)

    conserved_traces = get_conserved_traces(
        al_X, al_Y, outfile=os.path.join(params.outdir, "conserved_regions.csv")
    )

    # n_colors = 1 + len(Xs)
    # colors = px.colors.sample_colorscale("plotly3", [n / (n_colors - 1) for n in range(n_colors)])

    traces = get_traces(
        al_X,
        al_Y,
        edlib_X,
        edlib_Y,
        gt_X,
        gt_Y,
        Xs,
        Ys,
        bridges_X,
        bridges_Y,
        first_genes,
        second_genes,
        conserved_traces,
    )

    layout = go.Layout(
        title="Dotplot", xaxis_title=params.asm1_name, yaxis_title=params.asm2_name
    )
    fig = go.Figure(traces, layout)
    fig.update_xaxes(scaleanchor="x", scaleratio=1)
    fig.update_yaxes(range=[0, al_Y[-1] * 1.10], constrain="domain")
    fig.write_html(os.path.join(params.outdir, "dotplot.html"))


if __name__ == "__main__":
    main()
