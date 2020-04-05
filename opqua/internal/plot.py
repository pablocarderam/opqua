
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


def populationPlot(file_name, data, compartment='Infected', hosts=True, vectors=False,
                   num_top_populations=7, track_specific_populations=[], save_data_to_file="",
                   x_label='Time', y_label='Hosts', legend_title='Population', figsize=(8, 4), dpi=200, palette=cb_palette, stacked=False):
    '''Plot infection composition'''

    pops = populationsDf(data, compartment=compartment, hosts=hosts, vectors=vectors, num_top_populations=num_top_populations,
        track_specific_populations=track_specific_populations, save_to_file=save_data_to_file)

    if pops.shape[1] > 0:
        plt.figure(figsize=figsize, dpi=dpi) # make new figure
        ax = plt.subplot(1, 1, 1) # get axis

        if stacked:
            ax.stackplot( pops['Time'], pops.drop(columns='Time').transpose(), labels=pops.drop(columns='Time').columns, colors=palette )
        else:
            for i,c in enumerate(pops.columns[1:]):
                ax.plot( pops['Time'], pops[c], label=c, color=palette[i] )


        plt.xlabel(x_label) # labels
        plt.ylabel(y_label)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        handles, labels = ax.get_legend_handles_labels() # get legend
        plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), title=legend_title) # show it

        plt.savefig(file_name, bbox_inches='tight') # save
    else:
        ax = None
        print('Nothing to plot! Check your data.')

    return ax


def compartmentPlot(file_name, data, populations=[], hosts=True, vectors=False,
                    save_data_to_file="", x_label='Time', y_label='Hosts', legend_title='Compartment',
                    figsize=(8, 4), dpi=200, palette=cb_palette, stacked=False):
    '''Plot infection composition'''

    comp = compartmentDf(data, populations=populations, hosts=hosts,
        vectors=vectors, save_to_file=save_data_to_file)

    if comp.shape[1] > 0:
        plt.figure(figsize=figsize, dpi=dpi) # make new figure
        ax = plt.subplot(1, 1, 1) # get axis

        if stacked:
            ax.stackplot( comp['Time'], comp.drop(columns='Time').transpose(), labels=comp.drop(columns='Time').columns, colors=palette )
        else:
            for i,c in enumerate(comp.columns[1:]):
                ax.plot( comp['Time'], comp[c], label=c, color=palette[i] )


        plt.xlabel(x_label) # labels
        plt.ylabel(y_label)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        handles, labels = ax.get_legend_handles_labels() # get legend
        plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), title=legend_title) # show it

        plt.savefig(file_name, bbox_inches='tight') # save
    else:
        ax = None
        print('Nothing to plot! Check your data.')

    return ax


def compositionPlot(file_name, data, populations=[], organism='Pathogens', hosts=True, vectors=False,
                    num_top_genomes=7, track_specific_genomes=[], save_data_to_file="",
                    x_label='Time', y_label='Infections', legend_title='Genotype', figsize=(8, 4), dpi=200, palette=cb_palette, stacked=True):
    '''Plot infection composition'''

    comp = compositionDf(data, populations=populations, organism=organism, hosts=hosts,
        vectors=vectors, num_top_genomes=num_top_genomes, track_specific_genomes=track_specific_genomes,
        save_to_file=save_data_to_file)

    if comp.shape[1] > 0:
        plt.figure(figsize=figsize, dpi=dpi) # make new figure
        ax = plt.subplot(1, 1, 1) # get axis

        if stacked:
            ax.stackplot( comp['Time'], comp.drop(columns='Time').transpose(), labels=comp.drop(columns='Time').columns, colors=palette )
        else:
            for i,c in enumerate(comp.columns[1:]):
                ax.plot( comp['Time'], comp[c], label=c, color=palette[i] )


        plt.xlabel(x_label) # labels
        plt.ylabel(y_label)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        handles, labels = ax.get_legend_handles_labels() # get legend
        plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), title=legend_title) # show it

        plt.savefig(file_name, bbox_inches='tight') # save
    else:
        ax = None
        print('Nothing to plot! Check your data.')

    return ax
