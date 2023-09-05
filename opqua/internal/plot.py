
"""Contains graphmaking methods."""

### Imports ###
import copy as cp
import numpy as np # handle arrays
import pandas as pd # data wrangling
import matplotlib.pyplot as plt # plots
import seaborn as sns # pretty plots
import scipy.cluster.hierarchy as sp_hie
import scipy.spatial as sp_spa

from opqua.internal.data import saveToDf, populationsDf, compartmentDf, \
    compositionDf, pathogenDistanceDf

CB_PALETTE = ["#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]
    # www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    # http://jfly.iam.u-tokyo.ac.jp/color/
DEF_CMAP = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True, reverse=True)


def populationsPlot(
        file_name, data, compartment='Infected', hosts=True, vectors=False,
        num_top_populations=7, track_specific_populations=[],
        save_data_to_file="", x_label='Time', y_label='Infected hosts',
        legend_title='Population', legend_values=[], figsize=(8, 4), dpi=200,
        palette=CB_PALETTE, stacked=False):
    """Create plot with aggregated totals per population across time.

    Creates a line or stacked line plot with dynamics of a compartment
    across populations in the model, with one line for each population.

    Arguments:
        file_name (String): file path, name, and extension to save plot under.
        data (pandas DataFrame): dataframe with model history as produced by `saveToDf` function.

    Keyword arguments:
        compartment (String): subset of hosts/vectors to count totals of, can be either
            'Naive','Infected','Recovered', or 'Dead'. (default 'Infected')
        hosts (Boolean): whether to count hosts. Defaults to True.
        vectors (Boolean): whether to count vectors. Defaults to False.
        num_top_populations (int): how many populations to count separately and include
            as columns, remainder will be counted under column "Other"; if <0,
            includes all populations in model. Defaults to 7.
        track_specific_populations (list of Strings): contains IDs of specific populations to have
            as a separate column if not part of the top num_top_populations
            populations. Defaults to [].
        save_data_to_file (String): file path and name to save model plot data under, no
            saving occurs if empty string. Defaults to "".
        x_label(String): X axis title. Defaults to 'Time'.
        y_label (String): Y axis title. Defaults to 'Hosts'.
        legend_title (String): legend title. Defaults to 'Population'.
        legend_values (list of Strings): labels for each trace, if empty list, uses population IDs.
            Defaults to [].
        figsize (array-like of two ints): dimensions of figure. Defaults to (8,4).
        dpi (int): figure resolution. Defaults to 200.
        palette (list of color Strings): color palette to use for traces. Defaults to `CB_PALETTE`.
        stacked (Boolean): whether to draw a regular line plot instead of a stacked one. Defaults to False.

    Returns:
        axis object for plot with model population dynamics as described above.
    """

    pops = populationsDf(
        data, compartment=compartment, hosts=hosts, vectors=vectors,
        num_top_populations=num_top_populations,
        track_specific_populations=track_specific_populations,
        save_to_file=save_data_to_file
        )

    if pops.shape[1] > 0:
        plt.figure(figsize=figsize, dpi=dpi)
        ax = plt.subplot(1, 1, 1)

        if stacked:
            if len(legend_values) > 0:
                labs = legend_values
            else:
                labs = pops.drop(columns='Time').columns

            ax.stackplot(
                pops['Time'], pops.drop(columns='Time').transpose(),
                labels=labs, colors=palette
                )
        else:
            if len(legend_values) > 0:
                labs = legend_values
            else:
                labs = pops.drop(columns='Time').columns

            for i,c in enumerate(pops.columns[1:]):
                ax.plot(
                    pops['Time'], pops[c], label=labs[i], color=palette[i]
                    )

        plt.xlabel(x_label)
        plt.ylabel(y_label)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        handles, labels = ax.get_legend_handles_labels()
        plt.legend(
            handles, labels, loc='center left', bbox_to_anchor=(1, 0.5),
            title=legend_title
            )

        plt.savefig(file_name, bbox_inches='tight')
    else:
        ax = None
        print('Nothing to plot! Check your data.')

    return ax


def compartmentPlot(
        file_name, data, populations=[], hosts=True, vectors=False,
        save_data_to_file="", x_label='Time', y_label='Hosts',
        legend_title='Compartment', legend_values=[], figsize=(8, 4), dpi=200,
        palette=CB_PALETTE, stacked=False):
    """Create plot with num. of naive, susc., inf., rec. hosts/vectors vs. time.

    Creates a line or stacked line plot with dynamics of all compartments
    (naive, infected, recovered, dead) across selected populations in the model,
    with one line for each compartment.

    Arguments:
        file_name (String): file path, name, and extension to save plot under.
        data (pandas DataFrame): dataframe with model history as produced by `saveToDf` function.

    Keyword arguments:
        populations (list of Strings): IDs of populations to include in analysis; if empty, uses all
            populations in model. Defaults to [].
        hosts (Boolean): whether to count hosts. Defaults to True.
        vectors (Boolean): whether to count vectors. Defaults to False.
        save_data_to_file (String): file path and name to save model data under, no saving
            occurs if empty string. Defaults to "".
        x_label (String): X axis title. Defaults to 'Time'.
        y_label (String): Y axis title. Defaults to 'Hosts'.
        legend_title (String): legend title. Defaults to 'Population'.
        legend_values (list of Strings): labels for each trace, if empty list, uses population IDs. 
            Defaults to [].
        figsize (array-like of two ints): dimensions of figure. Defaults to (8,4).
        dpi (int): figure resolution. Defaults to 200.
        palette (list of color Strings): color palette to use for traces. Defaults to `CB_PALETTE`.
        stacked (Boolean): whether to draw a regular line plot instead of a stacked one. 
            Defaults to False.

    Returns:
        axis object for plot with model compartment dynamics as described above.
    """

    comp = compartmentDf(data, populations=populations, hosts=hosts,
        vectors=vectors, save_to_file=save_data_to_file)

    if comp.shape[1] > 0:
        plt.figure(figsize=figsize, dpi=dpi)
        ax = plt.subplot(1, 1, 1)

        if stacked:
            if len(legend_values) > 0:
                labs = legend_values
            else:
                labs = comp.drop(columns='Time').columns

            ax.stackplot(
                comp['Time'], comp.drop(columns='Time').transpose(),
                labels=labs, colors=palette
                )
        else:
            if len(legend_values) > 0:
                labs = legend_values
            else:
                labs = comp.drop(columns='Time').columns

            for i,c in enumerate(comp.columns[1:]):
                ax.plot(
                    comp['Time'], comp[c], label=labs[i], color=palette[i]
                    )

        plt.xlabel(x_label)
        plt.ylabel(y_label)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        handles, labels = ax.get_legend_handles_labels()
        plt.legend(
            handles, labels, loc='center left', bbox_to_anchor=(1, 0.5),
            title=legend_title
            )

        plt.savefig(file_name, bbox_inches='tight')
    else:
        ax = None
        print('Nothing to plot! Check your data.')

    return ax

def compositionPlot(
        file_name, data, composition_dataframe=None,
        populations=[], type_of_composition='Pathogens',
        hosts=True, vectors=False, num_top_sequences=7,
        track_specific_sequences=[], genomic_positions=[],
        count_individuals_based_on_model=None, save_data_to_file="",
        x_label='Time', y_label='Infections', legend_title='Genotype',
        legend_values=[], figsize=(8, 4), dpi=200, palette=CB_PALETTE,
        stacked=True, remove_legend=False, population_fraction=False, **kwargs):
    """Create plot with counts for pathogen genomes or resistance across time.

    Creates a line or stacked line plot with dynamics of the pathogen strains or
    protection sequences across selected populations in the model,
    with one line for each pathogen genome or protection sequence being shown.

    Of note: sum of totals for all sequences in one time point does not
    necessarily equal the number of infected hosts and/or vectors, given
    multiple infections in the same host/vector are counted separately.

    Arguments:
        file_name (String): file path, name, and extension to save plot under.
        data (pandas DataFrame): dataframe with model history as produced by `saveToDf` function.

    Keyword arguments:
        composition_dataframe (Pandas DataFrame): output of compositionDf() if already computed. 
            Defaults to None.
        populations (list of Strings): IDs of populations to include in analysis; if empty, uses all
            populations in model. Defaults to [].
        type_of_composition (String): field of data to count totals of, can be either
            'Pathogens' or 'Protection'. Defaults to 'Pathogens'.
        hosts (Boolean): whether to count hosts. Defaults to True.
        vectors (Boolean): whether to count vectors. Defaults to False.
        num_top_sequences (int): how many sequences to count separately and include
            as columns, remainder will be counted under column "Other"; if <0,
            includes all genomes in model. Defaults to 7.
        track_specific_sequences (list of Strings): contains specific sequences to have
            as a separate column if not part of the top num_top_sequences
            sequences. Defaults to [].
        genomic_positions (list of lists of int): list in which each element is a list with loci
            positions to extract (e.g. genomic_positions=[ [0,3], [5,6] ] extracts
            positions 0, 1, 2, and 5 from each genome); if empty, takes full genomes. Defaults to [].
        count_individuals_based_on_model (None or Model): Model object with populations and
            fitness functions used to evaluate the most fit pathogen genome in each
            host/vector in order to count only a single pathogen per host/vector, as
            opposed to all pathogens within each host/vector; if None, counts all
            pathogens. Defaults to None.
        save_data_to_file (String): file path and name to save model data under, no saving
            occurs if empty string. Defaults to "".
        x_label (String): X axis title. Defaults to 'Time'.
        y_label (String): Y axis title. Defaults to 'Hosts'.
        legend_title (String): legend title. Defaults to 'Population'.
        legend_values (list of Strings): labels for each trace, if empty list, uses population IDs. 
            Defaults to [].
        figsize (int): dimensions of figure. Defaults to (8,4).
        dpi (int): figure resolution. Defaults to 200.
        palette (list of color Strings): color palette to use for traces. Defaults to `CB_PALETTE`.
        stacked (Boolean): whether to draw a regular line plot instead of a stacked one. Defaults to False.
        remove_legend (Boolean): whether to print the sequences on the figure legend instead
            of printing them on a separate csv file. Defaults to True.
        population_fraction (Boolean): whether to graph fractions of pathogen population
            instead of pathogen counts. Defaults to False.
        **kwargs: additional arguents for joblib multiprocessing.

    Returns:
        axis object for plot with model sequence composition dynamics as described.
    """

    if composition_dataframe is None:
        comp = compositionDf(
            data, populations=populations,
            type_of_composition=type_of_composition,
            hosts=hosts, vectors=vectors, num_top_sequences=num_top_sequences,
            track_specific_sequences=track_specific_sequences,
            genomic_positions=genomic_positions,
            count_individuals_based_on_model=count_individuals_based_on_model,
            save_to_file=save_data_to_file, **kwargs
            )
    else:
        comp = composition_dataframe

    if population_fraction:
        total_infections = comp.drop(columns=['Time']).sum(axis=1)
        for col in comp.columns[1:]:
            print('corrected '+col)
            comp[col] = comp[col] / total_infections

    if comp.shape[1] > 1:
        plt.figure(figsize=figsize, dpi=dpi)
        ax = plt.subplot(1, 1, 1)

        if stacked:
            if len(legend_values) > 0:
                labs = legend_values
            else:
                labs = comp.drop(columns='Time').columns

            ax.stackplot(
                comp['Time'], comp.drop(columns='Time').transpose(),
                labels=labs, colors=palette
                )
        else:
            if len(legend_values) > 0:
                labs = legend_values
            else:
                labs = comp.drop(columns='Time').columns

            for i,c in enumerate(comp.columns[1:]):
                ax.plot(
                    comp['Time'], comp[c], label=labs[i], color=palette[i]
                    )

        plt.xlabel(x_label)
        plt.ylabel(y_label)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        if remove_legend:
            pd.DataFrame(
                labs, columns=['Groups']
                ).to_csv(file_name.split('.')[0]+'_labels.csv')
        else:
            handles, labels = ax.get_legend_handles_labels()
            plt.legend(
                handles, labels, loc='center left', bbox_to_anchor=(1, 0.5),
                title=legend_title
                )

        plt.savefig(file_name, bbox_inches='tight')
    else:
        ax = None
        print('Nothing to plot! Check your data.')

    return ax

def clustermap(
        file_name, data, num_top_sequences=-1, track_specific_sequences=[],
        seq_names=[], n_cores=0, method='weighted', metric='euclidean',
        save_data_to_file="", legend_title='Distance', legend_values=[],
        figsize=(10,10), dpi=200, color_map=DEF_CMAP):
    """Create a heatmap and dendrogram for pathogen genomes in data passed.

    Arguments:
        file_name (String): file path, name, and extension to save plot under.
        data (pandas DataFrame): dataframe with model history as produced by `saveToDf` function.

    Keyword arguments:
        num_top_sequences (int): how many sequences to include in matrix; if <0,
            includes all genomes in data passed. Deafults to -1.
        track_specific_sequences (list of Strings): contains specific sequences to include in matrix
            if not part of the top num_top_sequences sequences. Defaults to [].
        seq_names (list of Strings): list with names to be used for sequence labels in matrix must
            be of same length as number of sequences to be displayed; if empty,
            uses sequences themselves. Defaults to [].
        n_cores (int): number of cores to parallelize distance compute across, if 0, all
            cores available are used. Defaults to 0.
        method (String): clustering algorithm to use with seaborn clustermap. Defaults to 'weighted'.
        metric (String): distance metric to use with seaborn clustermap. Defaults to 'euclidean'.
        save_data_to_file (String): file path and name to save model data under, no saving
            occurs if empty string. Defaults to "".
        legend_title (String): legend title. Defaults to 'Distance'.
        figsize (array-like of two ints): dimensions of figure. Defaults to (8,4).
        dpi (int): figure resolution. Defaults to 200.
        color_map (cmap object): color map to use for traces. Defaults to `DEF_CMAP`.

    Returns:
        figure object for plot with heatmap and dendrogram as described.
    """

    dis = pathogenDistanceDf(
        data, num_top_sequences=num_top_sequences,
        track_specific_sequences=track_specific_sequences, seq_names=seq_names,
        n_cores=n_cores, save_to_file=save_data_to_file
        )

    lin = sp_hie.linkage(
        sp_spa.distance.squareform(dis), method='weighted',
        optimal_ordering=True
        )

    g = sns.clustermap(
        dis, method=method, metric=metric, cbar_kws={'label': legend_title},
        cmap=color_map, figsize=figsize
        )

    g.savefig(file_name, dpi=dpi, bbox_inches='tight')

    return g
