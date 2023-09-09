import pandas as pd
import numpy as np


def projection(g, al_S, al_D, al_i):
    if al_S[al_i] == al_S[al_i - 1]:
        return al_D[al_i - 1]
    return al_D[al_i - 1] + (g - al_S[al_i - 1]) * (al_D[al_i] - al_D[al_i - 1]) / (
        al_S[al_i] - al_S[al_i - 1]
    )


def get_projections(genes_starts_df, al_S, al_D):
    orthologue_pos = []
    # iterate over alignment to find orthologue pos
    al_i = 1  # index of upper bound of current interval
    genes_i = 0
    while genes_i < len(genes_starts_df):
        pos = genes_starts_df.pos[genes_i]
        if al_i == len(al_S):
            # if pos==al_S[al_i-1]: # gene at the end of alignment
            # # else gene falls out of alignment (not ideal, TODO add warning)
            orthologue_pos.append(al_D[-1])
            genes_i += 1
        elif pos < al_S[al_i]:
            # gene contained in this interval
            orthologue_pos.append(projection(pos, al_S, al_D, al_i))
            genes_i += 1
        else:
            # gene ahead of this interval
            al_i += 1
    return orthologue_pos


def get_first_genes(genes_starts_csv, al_X, al_Y):
    genes_starts_df = pd.read_csv(genes_starts_csv)
    genes_starts_df = genes_starts_df.sort_values("pos").reset_index()
    orthologue_pos = get_projections(genes_starts_df, al_X, al_Y)
    colors = ["gray"] * len(genes_starts_df)
    if "color" in genes_starts_df.columns:
        colors = genes_starts_df.color
    return [
        ([p, p], [0, o], [p, 0], [o, o], name, color)
        for p, o, name, color in zip(
            genes_starts_df.pos, orthologue_pos, genes_starts_df.gene, colors
        )
    ]


def get_second_genes(genes_starts_csv, al_X, al_Y):
    genes_starts_df = pd.read_csv(genes_starts_csv)
    genes_starts_df = genes_starts_df.sort_values("pos").reset_index()
    orthologue_pos = get_projections(genes_starts_df, al_Y, al_X)
    colors = ["maroon"] * len(genes_starts_df)
    if "color" in genes_starts_df.columns:
        colors = genes_starts_df.color
    return [
        ([o, 0], [p, p], [o, o], [0, p], name, color)
        for p, o, name, color in zip(
            genes_starts_df.pos, orthologue_pos, genes_starts_df.gene, colors
        )
    ]
