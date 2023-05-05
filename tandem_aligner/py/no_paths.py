# TODO: add labels

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    params = parser.parse_args()
    no_paths_file =params.input # "/home/i3gupta/tandem_aligner/tandem_aligner/test_dataset/fig1/result/no_paths.csv"

    df = pd.read_csv(no_paths_file)

    print(df)

    df["n_paths"] = df["n_opt_paths_fwd"]* df["n_opt_paths_rev"]
    df["pct_paths"] = (df["n_paths"]+10000)/(df["n_paths"].max()+10000)

    print("Total paths found = ", df["n_paths"].max())

    #TOREMOVE
    print(df["pct_paths"]) 

    fig,axs = plt.subplots(3)
    fig.set_size_inches(10,10*len(axs))
    axs[0].xaxis.tick_top()
    axs[1].xaxis.tick_top()
    axs[2].xaxis.tick_top()
    axs[0].invert_yaxis()
    axs[1].invert_yaxis()
    axs[2].invert_yaxis()

    axs[1].set_title("Bridges")

    for i,row in df.iterrows():
        # graph 1
        axs[0].plot([row.start_2,row.start_2 + row.len],[row.start_1,row.start_1+row.len], c = str(1-row.pct_paths), alpha=0.5)
        # graph 2
        if (row.pct_paths ==1):
            axs[1].plot([row.start_2,row.start_2 + row.len],[row.start_1,row.start_1+row.len], c = str(1-row.pct_paths))
        # graph 3
        if (row.pct_paths ==1):
            c = 'red'
        elif(row.n_paths == 0):
            c = 'green'
        elif(row.n_count==1 and row.m_count==1): 
            c = 'blue'
        else: 
            c = 'yellow'
        axs[2].plot([row.start_2,row.start_2 + row.len],[row.start_1,row.start_1+row.len], c = c, alpha = 0.5)
    fig.savefig("/".join(os.path.split(no_paths_file)[:-1]) + "/no_paths_plot")
    fig = plt.figure()
    plt.hist(df["n_paths"])
    fig.savefig("/".join(os.path.split(no_paths_file)[:-1]) + "/no_paths_hist")

main()