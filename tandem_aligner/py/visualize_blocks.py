import argparse
from collections import namedtuple
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from utils.os_utils import expandpath, smart_makedirs


AlignedBlock = namedtuple("AlignedBlock", ["st1", "en1", "st2", "en2"])


# taken from https://stackoverflow.com/a/32520273
def rand_cmap(
    nlabels, type="bright", first_color_black=True, last_color_black=False, verbose=True
):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    if type not in ("bright", "soft"):
        print('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print("Number of labels: " + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == "bright":
        randHSVcolors = [
            (
                np.random.uniform(low=0.0, high=1),
                np.random.uniform(low=0.2, high=1),
                np.random.uniform(low=0.9, high=1),
            )
            for i in range(nlabels)
        ]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(
                colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2])
            )

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list(
            "new_map", randRGBcolors, N=nlabels
        )

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == "soft":
        low = 0.6
        high = 0.95
        randRGBcolors = [
            (
                np.random.uniform(low=low, high=high),
                np.random.uniform(low=low, high=high),
                np.random.uniform(low=low, high=high),
            )
            for i in range(nlabels)
        ]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list(
            "new_map", randRGBcolors, N=nlabels
        )

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(
            ax,
            cmap=random_colormap,
            norm=norm,
            spacing="proportional",
            ticks=None,
            boundaries=bounds,
            format="%1i",
            orientation="horizontal",
        )

    return random_colormap, randRGBcolors


def read_blocks(fn):
    blocks = []
    with open(fn) as f:
        for line in f:
            st1, en1, st2, en2 = [int(x) for x in line.split()]
            blocks.append(AlignedBlock(st1, en1, st2, en2))
    return blocks


def plot_blocks(blocks, str1_name, str2_name, outdir, min_len=1):
    fig, ax = plt.subplots(2, figsize=(25, 8), dpi=300)

    # _, color = rand_cmap(len(blocks), type='bright', first_color_black=True, last_color_black=False, verbose=True)
    # print(color)

    color = plt.cm.rainbow(np.linspace(0, 1, len(blocks)))
    np.random.seed(0)
    np.random.shuffle(color)

    for block, col in zip(blocks, color):
        rect1 = matplotlib.patches.Rectangle(
            (block.st1, 0), block.en1 - block.st1, 1, color=col, alpha=0.7
        )
        rect2 = matplotlib.patches.Rectangle(
            (block.st2, 0), block.en2 - block.st2, 1, color=col, alpha=0.7
        )

        ax[0].add_patch(rect1)
        ax[1].add_patch(rect2)

    xmax = max(blocks[-1].en1, blocks[-1].en2)
    ax[0].set_xlim(0, xmax)
    ax[1].set_xlim(0, xmax)
    step = 10 ** (int(np.log10(xmax)) - 1)
    ax[1].set_xticks(np.arange(0, xmax, step))
    ax[1].set_xticklabels(labels=np.arange(0, xmax, step), rotation=30, fontsize=14)

    ax[0].tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False,
    )  # labels along the bottom edge are off

    ax[0].tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    ax[1].tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    ax[0].set_ylabel(str1_name[:10], fontsize=20)
    ax[1].set_ylabel(str2_name[:10], fontsize=20)
    outfn = os.path.join(outdir, "blocks.pdf")
    plt.savefig(outfn, format="pdf")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--blocks", required=True)
    parser.add_argument("-1", "--name1", required=True)
    parser.add_argument("-2", "--name2", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    params = parser.parse_args()

    params.outdir = expandpath(params.outdir)
    smart_makedirs(params.outdir)

    blocks = read_blocks(params.blocks)
    print(blocks)
    plot_blocks(blocks, params.name1, params.name2, params.outdir)


if __name__ == "__main__":
    main()
