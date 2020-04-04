# TODO: write data wrangling for time series, graphs, phylogeny
# TODO: distance matrix
# TODO: gráfica igual pero para inmunidad,
# TODO: gráficas "normales" de compartimientos SIRD dentro de una sola población
# TODO: gráficas de infectados por población
# TODO: filogenias de la infección

### Imports ###
import matplotlib.pyplot as plt # plots
import seaborn as sns # pretty plots
import numpy as np # handle arrays
import pandas as pd # data wrangling
import copy as cp

from opqua.internal.data import *

cb_palette = ["#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]
    # www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    # http://jfly.iam.u-tokyo.ac.jp/color/

def compositionPlot(file_name, data, populations=[], type='Pathogens', hosts=True, vectors=True,
                    num_top_genomes=7, track_specific_genomes=[], save_data_to_file="",
                    x_label='Time', y_label='Infections', figsize=(8, 4), dpi=200, palette=cb_palette):
    '''Plot infection composition'''

    comp = compositionDf(data, populations=populations, type=type, hosts=hosts,
        vectors=vectors, num_top_genomes=num_top_genomes, track_specific_genomes=track_specific_genomes,
        save_to_file=save_data_to_file)

    if comp.shape[1] > 0:
        plt.figure(figsize=figsize, dpi=dpi) # make new figure
        ax = plt.subplot(1, 1, 1) # get axis

        ax.stackplot(comp.index, comp.transpose(), labels=comp.columns, colors=palette)

        plt.xlabel(x_label) # labels
        plt.ylabel(y_label)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        handles, labels = ax.get_legend_handles_labels() # get legend
        plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5)) # show it

        plt.savefig(file_name, bbox_inches='tight') # save
    else:
        ax = None
        print('Nothing to plot! Check your data.')

    return ax