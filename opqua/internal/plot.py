# TODO: write data wrangling for time series, graphs, phylogeny
# TODO: distance matrix

### Imports ###
import matplotlib.pyplot as plt # plots
import seaborn as sns # pretty plots
import numpy as np # handle arrays
import pandas as pd # data wrangling

cb_palette = ["#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]
    # www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    # http://jfly.iam.u-tokyo.ac.jp/color/

def infectionsComposition(file_name, data, populations=[], hosts=True,
                          num_top_genomes=8, track_specific_genotypes=[],
                          x_label='Time', y_label='Infections',
                          figsize=(8, 4), dpi=200, palette=cb_palette):
    '''Plot infection composition'''

    if len(populations) > 0:
        dat = data[ data['Population'].isin( populations ) ]
    else:
        dat = data

    all_genomes = pd.Series( ';'.join( dat['Pathogens'].dropna() ).split(';') ).str.strip()
    top_genomes = all_genomes.value_counts(ascending=False)
    if len(top_genomes) > num_top_genomes:
        top_genomes = top_genomes[0:num_top_genomes]

    genome_split_data = []
    for genome in top_genomes.index:
        dat_genome = dat[ dat['Pathogens'].str.contains(genome, na=False) ]
        grouped = dat_genome.groupby('Time').size().reset_index(name=genome)
        grouped = grouped.set_index('Time')
        genome_split_data.append(grouped)

    times = pd.DataFrame( index=pd.unique( data['Time'] ) )
    comp = times.join( genome_split_data, how='outer' )

    plt.figure(figsize=figsize, dpi=dpi) # make new figure
    ax = plt.subplot(1, 1, 1) # get axis

    ax.stackplot(comp.index, comp.transpose(), labels=top_genomes.index, colors=palette)

    plt.xlabel(x_label) # labels
    plt.ylabel(y_label)
    handles, labels = ax.get_legend_handles_labels() # get legend
    plt.legend(handles, labels, loc='upper right') # show it

    plt.savefig('file_name', bbox_inches='tight') # save

    return ax
