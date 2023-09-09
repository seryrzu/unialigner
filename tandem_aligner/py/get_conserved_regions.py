import numpy as np


def get_conserved_regions(al_X, al_Y, max_indel_size, min_conserved_size):
    start_X, start_Y = 0, 0
    conserved_regions = []
    for i in range(1, len(al_X)):
        if (al_X[i] == al_X[i - 1] or al_Y[i] == al_Y[i - 1]) and (
            max(al_X[i] - al_X[i - 1], al_Y[i] - al_Y[i - 1]) > max_indel_size
        ):
            if min(al_X[i - 1] - start_X, al_Y[i - 1] - start_Y) > min_conserved_size:
                conserved_regions.append(
                    ((start_X, start_Y), (al_X[i - 1], al_Y[i - 1]))
                )
            start_X, start_Y = al_X[i], al_Y[i]
    return conserved_regions


def write_conserved_regions(conserved_regions, outfile):
    out = open(outfile, "w")
    out.write("first_start,second_start,first_end,second_end\n")
    for conserved_region in conserved_regions:
        out.write(
            ",".join([str(x) for x in conserved_region[0] + conserved_region[1]]) + "\n"
        )


def get_conserved_traces(
    al_X, al_Y, max_indel_size=1e3, min_conserved_size=1e3, outfile=None
):
    conserved_regions = get_conserved_regions(
        al_X, al_Y, max_indel_size, min_conserved_size
    )
    if outfile is not None:
        write_conserved_regions(conserved_regions, outfile)
    conserved_traces = []
    for i in range(len(conserved_regions)):
        conserved_traces.extend(conserved_regions[i])
        conserved_traces.extend([(np.nan, np.nan)])
    return np.array(conserved_traces)
