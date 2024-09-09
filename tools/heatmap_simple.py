#! /data/cqs/chenh19/project/nrsa_v2/miniconda3/bin/python3.12
import pandas as pd
import numpy as np
import os, sys, re

def draw_heatmap_pindex(pwout, fdr_thres=0.05, fc_thres=0):
    """
    plot for the pindex change
    fc_thres, the fc absolution value threshold for the pindex change
    """
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap
    
    fn_pindex_change = f'{pwout}/known_gene/pindex_change.txt'
    fno = f'pindex_change.pdf'
    if not os.path.exists(fn_pindex_change):
        logger.error(f'{fn_pindex_change} not found')
        return 1

    data_plot = pd.read_csv(fn_pindex_change, sep='\t')
    # Transcript	Gene	log2fc	pvalue	FDR
    col_fdr, col_logfc = 'FDR', 'log2fc'
    data_plot[col_logfc] = data_plot[col_logfc].replace([np.inf, -np.inf], np.nan)

    # drop the rows which have NA in the log2fc and fdr column
    data_plot = data_plot.dropna(subset=[col_logfc, col_fdr])

    # filter out the genes with FDR > 0.05
    data_plot = data_plot[(data_plot[col_fdr] < fdr_thres) & (abs(data_plot[col_logfc]) > fc_thres)]
    # sort by log2FC
    data_plot = data_plot.sort_values(by=col_logfc)

    with PdfPages(fno) as pdf:
        plt.figure(figsize=(2, 6))
        grid = plt.GridSpec(2, 1, height_ratios=[14, 1])
        ax1 = plt.subplot(grid[0, 0])
        values = data_plot['log2fc'].values[::-1]
        cmap = LinearSegmentedColormap.from_list("my_palette", ["green", "yellow", "red"], N=209)
        ax1.imshow(values.reshape(-1, 1), cmap=cmap, aspect='auto')
        ax1.axis('off')
        
        
        ax2 = plt.subplot(grid[1, 0])
        gradient = np.linspace(min(values), max(values), 209).reshape(1, -1)
        ax2.imshow(gradient, aspect='auto', cmap=cmap)
        ax2.set_xticks([0, 104, 208])
        ax2.set_xticklabels([round(min(values), 2), "0", round(max(values), 2)], fontsize=8)
        ax2.set_yticks([])
        ax2.set_xlabel("")
        ax2.set_ylabel("")
        plt.subplots_adjust(left = 0, right = 0.99, top=0.95, bottom = 0, hspace=0.05)
        pdf.savefig(bbox_inches='tight')
        plt.close()


import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('pwout', help="""output folder with known_gene and intermediate""")
args = ps.parse_args()

draw_heatmap_pindex(args.pwout)

