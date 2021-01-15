
"""Contains data wrangling methods."""

import numpy as np # handle arrays
import pandas as pd # data wrangling
import copy as cp
import joblib as jl
import textdistance as td
import scipy.spatial.distance as sp_dist

def saveToDf(history,save_to_file,n_cores=0):
    """Save status of model to dataframe, write to file location given.

    Creates a pandas Dataframe in long format with the given model history, with
    one host or vector per simulation time in each row, and columns:
        Time - simulation time of entry
        Population - ID of this host/vector's population
        Organism - host/vector
        ID - ID of host/vector
        Pathogens - all genomes present in this host/vector separated by ;
        Protection - all genomes present in this host/vector separated by ;
        Alive - whether host/vector is alive at this time, True/False

    Writing straight to a file and then reading into a pandas dataframe was
    actually more efficient than concatenating directly into a pd dataframe.

    Arguments:
    history -- dictionary containing model state history, with keys=times and
        values=Model objects with model snapshot at that time point
    save_to_file -- file path and name to save model data under (String)

    Keyword arguments:
    n_cores -- number of cores to parallelize file export across, if 0, all
        cores available are used (default 0; int)

    Returns:
    pandas dataframe with model history as described above
    """

    print('Saving file...')

    if not n_cores:
        n_cores = jl.cpu_count()

    new_df = ','.join(
        ['Time','Population','Organism','ID','Pathogens','Protection','Alive']
        ) + '\n' + '\n'.join( jl.Parallel(n_jobs=n_cores, verbose=10) (
            jl.delayed( lambda d: ''.join(d) ) (
                '\n'.join( [
                    '\n'.join( [ ','.join( [
                        str(time), str(pop.id), 'Host', str(host.id), '"'
                        + ';'.join( host.pathogens.keys() )
                        + '"', '"' + ';'.join( host.protection_sequences )
                        + '"', 'True'
                        ] ) for host in pop.hosts ] ) + '\n'
                    + '\n'.join( [ ','.join( [
                        str(time), str(pop.id), 'Vector', str(vector.id), '"'
                        + ';'.join( vector.pathogens.keys() )
                        + '"', '"' + ';'.join( vector.protection_sequences )
                        + '"', 'True'
                        ] ) for vector in pop.vectors ] ) + '\n'
                    + '\n'.join( [ ','.join( [
                        str(time), str(pop.id), 'Host', str(host.id), '"'
                        + ';'.join( host.pathogens.keys() ) + '"', '"'
                        + ';'.join( host.protection_sequences ) + '"', 'False'
                        ] ) for host in pop.dead_hosts ] ) + '\n'
                    + '\n'.join( [ ','.join( [
                        str(time), str(pop.id), 'Vector', str(vector.id), '"'
                        + ';'.join( vector.pathogens.keys() ) + '"', '"'
                        + ';'.join( vector.protection_sequences ) + '"', 'False'
                        ] ) for vector in pop.dead_vectors ] )
                    for id,pop in model.populations.items()
                ] )
            ) for time,model in history.items()
        ) )

    new_df = new_df.replace(
        '\n\n','\n'
        ).replace('\n\n','\n').replace('\n\n','\n')

    file = open(save_to_file,'w')
    file.write(new_df)
    file.close()

    new_df = pd.read_csv(save_to_file)

    print('...file saved.')

    return new_df

def populationsDf(
        data, compartment='Infected', hosts=True, vectors=False,
        num_top_populations=-1, track_specific_populations=[], save_to_file=""):
    """Create dataframe with aggregated totals per population.

    Creates a pandas Dataframe in long format with dynamics of a compartment
    across populations in the model, with one time point in each row and columns
    for time as well as each population.

    Arguments:
    data -- dataframe with model history as produced by saveToDf function

    Keyword arguments:
    compartment -- subset of hosts/vectors to count totals of, can be either
        'Naive','Infected','Recovered', or 'Dead' (default 'Infected'; String)
    hosts -- whether to count hosts (default True, Boolean)
    vectors -- whether to count vectors (default False, Boolean)
    num_top_populations -- how many populations to count separately and include
        as columns, remainder will be counted under column "Other"; if <0,
        includes all populations in model (default -1; int)
    track_specific_populations -- contains IDs of specific populations to have
        as a separate column if not part of the top num_top_populations
        populations (list of Strings)
    save_to_file -- file path and name to save model data under, no saving
        occurs if empty string (default ''; String)

    Returns:
    pandas dataframe with model population dynamics as described above
    """

    dat = cp.deepcopy( data )

    if not hosts:
        dat = dat[ dat['Organism'] != 'Host' ]
    if not vectors:
        dat = dat[ dat['Organism'] != 'Vector' ]

    if num_top_populations < 0:
        num_top_populations = len( pd.unique( dat['Population'] ) )

    dat['Infected'] = ( dat['Pathogens'].fillna('').str.len() > 0 )
    dat['Protected'] = ( dat['Protection'].fillna('').str.len() > 0 )

    grouped = dat.groupby( [
        'Time','Population','Alive','Infected','Protected'
        ] ).size().reset_index(name='Number')

    compartment_names = ['Naive','Infected','Recovered','Dead']

    grouped['Compartment'] = compartment_names[3]
    grouped.loc[
        ( grouped['Alive'] == True ) & ( grouped['Infected'] == False )
            & ( grouped['Protected'] == False ), 'Compartment'
        ] = compartment_names[0]
    grouped.loc[
        ( grouped['Alive'] == True ) & ( grouped['Infected'] == True ),
        'Compartment'
        ] = compartment_names[1]
    grouped.loc[
        ( grouped['Alive'] == True ) & ( grouped['Infected'] == False )
        & ( grouped['Protected'] == True ), 'Compartment'
        ] = compartment_names[2]

    grouped = grouped[ grouped['Compartment'] == compartment ]

    grouped = grouped.groupby(
        ['Time','Population','Compartment']
    ).sum().reset_index()

    grouped = grouped.drop(
        columns=['Alive','Infected','Protected','Compartment']
        )

    grouped = grouped.pivot(
        columns='Population', values='Number', index='Time'
        ).fillna(0).reset_index('Time')

    for pop in pd.unique(data['Population']):
        if pop not in grouped.columns:
            grouped[pop] = 0

    if len(grouped.columns)-1 < num_top_populations:
        num_top_populations = len(grouped.columns)-1

    populations_to_drop = list(grouped.columns)[ num_top_populations+1: ]
    for pop in track_specific_populations:
        if pop in populations_to_drop:
            populations_to_drop.remove(pop)

    grouped['Other'] = 0
    if len(populations_to_drop) > 0:
        grouped['Other'] = grouped[populations_to_drop].sum(axis=1)

    if grouped['Other'].sum() == 0:
        populations_to_drop = populations_to_drop + ['Other']

    grouped = grouped.drop( columns=populations_to_drop )

    if len(save_to_file) > 0:
        grouped.to_csv(save_to_file, index=False)

    return grouped


def compartmentDf(
        data, populations=[], hosts=True, vectors=False, save_to_file=""):
    """Create dataframe with number of naive, susc., inf., rec. hosts/vectors.

    Creates a pandas Dataframe with dynamics of all compartments (naive,
    infected, recovered, dead) across selected populations in the model,
    with one time point in each row and columns for time as well as each
    compartment.

    Arguments:
    data -- dataframe with model history as produced by saveToDf function

    Keyword arguments:
    populations -- IDs of populations to include in analysis; if empty, uses all
        populations in model (default empty list; list of Strings)
    hosts -- whether to count hosts (default True, Boolean)
    vectors -- whether to count vectors (default False, Boolean)
    save_to_file -- file path and name to save model data under, no saving
        occurs if empty string (default ''; String)

    Returns:
    pandas dataframe with model compartment dynamics as described above
    """

    if len(populations) > 0:
        dat = cp.deepcopy( data[ data['Population'].isin( populations ) ] )
    else:
        dat = cp.deepcopy( data )

    if not hosts:
        dat = dat[ dat['Organism'] != 'Host' ]

    if not vectors:
        dat = dat[ dat['Organism'] != 'Vector' ]

    dat['Infected'] = ( dat['Pathogens'].fillna('').str.len() > 0 )
    dat['Protected'] = ( dat['Protection'].fillna('').str.len() > 0 )

    grouped = dat.groupby( [
        'Time','Alive','Infected','Protected'
        ] ).size().reset_index(name='Number')

    compartment_names = ['Naive','Infected','Recovered','Dead']

    grouped['Compartment'] = compartment_names[3]
    grouped.loc[
        ( grouped['Alive'] == True ) & ( grouped['Infected'] == False )
        & ( grouped['Protected'] == False ), 'Compartment'
        ] = compartment_names[0]
    grouped.loc[
        ( grouped['Alive'] == True ) & ( grouped['Infected'] == True )
        & ( grouped['Protected'] == True ), 'Compartment'
        ] = compartment_names[1]
    grouped.loc[
        ( grouped['Alive'] == True ) & ( grouped['Infected'] == True )
        & ( grouped['Protected'] == False ), 'Compartment'
        ] = compartment_names[1] + '_2'
    grouped.loc[
        ( grouped['Alive'] == True ) & ( grouped['Infected'] == False )
        & ( grouped['Protected'] == True ), 'Compartment'
        ] = compartment_names[2]

    grouped = grouped.drop( columns=['Alive','Infected','Protected'] )
    grouped = grouped.pivot(
        columns='Compartment', values='Number', index='Time'
        ).fillna(0).reset_index('Time')

    if ( compartment_names[1] in grouped.columns
            and compartment_names[1] + '_2' in grouped.columns ):
        grouped[compartment_names[1]] = ( grouped[ compartment_names[1] ]
            + grouped[ compartment_names[1] + '_2' ] )
        grouped = grouped.drop( columns=[ compartment_names[1] + '_2' ])
    elif ( compartment_names[1] + '_2' in grouped.columns ):
        grouped[compartment_names[1]] = grouped[ compartment_names[1] + '_2' ]
        grouped = grouped.drop( columns=[ compartment_names[1] + '_2' ])

    for comp_name in compartment_names:
        if comp_name not in grouped.columns:
            grouped[comp_name] = 0


    if len(save_to_file) > 0:
        grouped.to_csv(save_to_file, index=False)

    return grouped


def compositionDf(
        data, populations=[], type_of_composition='Pathogens', hosts=True,
        vectors=False, num_top_sequences=-1, track_specific_sequences=[],
        save_to_file=""):
    """Create dataframe with counts for pathogen genomes or resistance.

    Creates a pandas Dataframe with dynamics of the pathogen strains or
    protection sequences across selected populations in the model,
    with one time point in each row and columns for pathogen genomes or
    protection sequences.

    Of note: sum of totals for all sequences in one time point does not
    necessarily equal the number of infected hosts and/or vectors, given
    multiple infections in the same host/vector are counted separately.

    Arguments:
    data -- dataframe with model history as produced by saveToDf function

    Keyword arguments:
    populations -- IDs of populations to include in analysis; if empty, uses all
        populations in model (default empty list; list of Strings)
    type_of_composition -- field of data to count totals of, can be either
        'Pathogens' or 'Protection' (default 'Pathogens'; String)
    hosts -- whether to count hosts (default True, Boolean)
    vectors -- whether to count vectors (default False, Boolean)
    num_top_sequences -- how many sequences to count separately and include
        as columns, remainder will be counted under column "Other"; if <0,
        includes all genomes in model (default -1; int)
    track_specific_sequences -- contains specific sequences to have
        as a separate column if not part of the top num_top_sequences
        sequences (default empty list; list of Strings)
    save_to_file -- file path and name to save model data under, no saving
        occurs if empty string (default ''; String)

    Returns:
    pandas dataframe with model sequence composition dynamics as described above
    """

    if len(populations) > 0:
        dat = cp.deepcopy( data[ data['Population'].isin( populations ) ] )
    else:
        dat = cp.deepcopy( data )

    if not hosts:
        dat = dat[ dat['Organism'] != 'Host' ]
    if not vectors:
        dat = dat[ dat['Organism'] != 'Vector' ]

    if num_top_sequences < 0:
        num_top_sequences = len( pd.unique( dat[type_of_composition] ) )

    all_genomes = pd.Series(
        ';'.join( dat[type_of_composition].dropna() ).split(';')
        ).str.strip()
    top_genomes = all_genomes.value_counts(ascending=False)

    if len(top_genomes) < num_top_sequences:
        num_top_sequences = len(top_genomes)

    genomes_to_track = list(top_genomes[0:num_top_sequences].index)
    for genome in track_specific_sequences:
        if genome not in genomes_to_track:
            genomes_to_track.append(genome)


    genome_split_data = []
    other_genome_data = pd.DataFrame(
        np.zeros( len( pd.unique( data['Time'] ) ) ),
        index=pd.unique( data['Time'] ), columns=['Other']
        )
    c = 0

    if len( ''.join(top_genomes.index) ) > 0:
        for genome in top_genomes.index:
            dat_genome = dat[ dat[type_of_composition].str.contains(
                genome, na=False
                ) ]
            grouped = dat_genome.groupby('Time').size().reset_index(name=genome)
            grouped = grouped.set_index('Time')

            if genome in genomes_to_track:
                genome_split_data.append(grouped)
            else:
                other_genome_data = other_genome_data.join(
                    grouped, how='outer'
                    ).fillna(0)
                other_genome_data['Other'] = ( other_genome_data['Other']
                    + other_genome_data[genome] )
                other_genome_data = other_genome_data.drop(columns=[genome])

            c += 1
            print(
                str(c) + ' / ' + str( len(top_genomes.index) )
                + ' genotypes processed.'
                )


        if other_genome_data['Other'].sum() > 0:
            genome_split_data += [other_genome_data]
            genomes_to_track += ['Other']


    times = pd.DataFrame(
        np.zeros( len( pd.unique( data['Time'] ) ) ),
        index=pd.unique( data['Time'] ), columns=['*None*']
        )
    composition = times.join( genome_split_data, how='outer' )
    composition = composition.drop(columns=['*None*']).fillna(0)
    composition = composition.reset_index()
    composition.columns = ['Time'] + list( composition.columns )[1:]

    if len(save_to_file) > 0:
        composition.to_csv(save_to_file, index=False)

    return composition


def getPathogens(data, save_to_file=""):
    """Create Dataframe with counts for all pathogen genomes in data.

    Returns sorted pandas Dataframe with counts for occurrences of all pathogen
    genomes in data passed.

    Arguments:
    data -- dataframe with model history as produced by saveToDf function

    Keyword arguments:
    save_to_file -- file path and name to save model data under, no saving
        occurs if empty string (default ''; String)

    Returns:
    pandas dataframe with Series as described above
    """

    out = pd.Series( ';'.join(
        data['Pathogens'].dropna()
        ).split(';') ).str.strip().value_counts(ascending=False).reset_index()
    out.columns = ['Pathogens','Counts']

    if len(save_to_file) > 0:
        out.to_csv(save_to_file, index=False)

    return out

def getProtections(data, save_to_file=""):
    """Create Dataframe with counts for all protection sequences in data.

    Returns sorted pandas Dataframe with counts for occurrences of all
    protection sequences in data passed.

    Arguments:
    data -- dataframe with model history as produced by saveToDf function

    Keyword arguments:
    save_to_file -- file path and name to save model data under, no saving
        occurs if empty string (default ''; String)

    Returns:
    pandas dataframe with Series as described above
    """

    out = pd.Series( ';'.join(
        data['Protection'].dropna()
        ).split(';') ).str.strip().value_counts(ascending=False).reset_index()
    out.columns = ['Protection','Counts']

    if len(save_to_file) > 0:
        out.to_csv(save_to_file, index=False)

    return out

def pathogenDistanceDf(
        data, num_top_sequences=-1, track_specific_sequences=[], seq_names=[],
        save_to_file="", n_cores=0):
    """Create DataFrame with pairwise distances for pathogen sequences in data.

    DataFrame has indexes and columns named according to genomes or argument
    seq_names, if passed. Distance is measured as percent Hamming distance from
    an optimal genome sequence.

    Arguments:
    data -- dataframe with model history as produced by saveToDf function

    Keyword arguments:
    num_top_sequences -- how many sequences to include in matrix; if <0,
        includes all genomes in data passed (default -1; int)
    track_specific_sequences -- contains specific sequences to include in matrix
        if not part of the top num_top_sequences sequences (default empty list;
        list of Strings)
    seq_names -- list with names to be used for sequence labels in matrix must
        be of same length as number of sequences to be displayed; if empty,
        uses sequences themselves (default empty list; list of Strings)
    save_to_file -- file path and name to save model data under, no saving
        occurs if empty string (default ''; String)
    n_cores -- number of cores to parallelize distance compute across, if 0, all
        cores available are used (default 0; int)

    Returns:
    pandas dataframe with distance matrix as described above
    """

    sequences = getPathogens(data)['Pathogens']

    if num_top_sequences > 0:
        sequences = sequences[0:num_top_sequences]

    str_mat = [ list(seq) for seq in sequences ]

    if not n_cores:
        n_cores = jl.cpu_count()


    dis_mat = np.array([[td.hamming(s1, s2) / len(s1)
        for s2 in sequences] for s1 in sequences])

    # For some reason, this code triggers a full rerun of the simulation.
    # Possible joblib bug?
    # dis_mat = np.array( jl.Parallel(n_jobs=n_cores, verbose=1) (
    #     jl.delayed( lambda s1: [
    #         td.hamming(s1, s2) / len(s1) for s2 in sequences
    #         ] ) (s1) for s1 in sequences
    #     ) )

    names = sequences
    if len(seq_names) > 0:
        names = seq_names

    dis_df = pd.DataFrame( dis_mat, index=names, columns=names )

    if len(save_to_file) > 0:
        dis_df.to_csv(save_to_file, index=True)

    return dis_df
