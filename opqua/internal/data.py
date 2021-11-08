
"""Contains data wrangling methods."""

import numpy as np # handle arrays
import pandas as pd # data wrangling
import copy as cp
import joblib as jl
import textdistance as td
import scipy.spatial.distance as sp_dist

def saveToDf(history,save_to_file,n_cores=0,verbose=10, **kwargs):
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
    **kwargs -- additional arguents for joblib multiprocessing

    Returns:
    pandas dataframe with model history as described above
    """

    print('Saving file...')

    if not n_cores:
        n_cores = jl.cpu_count()

    new_df = ','.join(
        ['Time','Population','Organism','ID','Pathogens','Protection','Alive']
        ) + '\n' + '\n'.join( jl.Parallel(
            n_jobs=n_cores, verbose=verbose, **kwargs) (
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
        genomic_positions=[], count_individuals_based_on_model=None,
        save_to_file="", n_cores=0, **kwargs):
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
    genomic_positions -- list in which each element is a list with loci
        positions to extract (e.g. genomic_positions=[ [0,3], [5,6] ] extracts
        positions 0, 1, 2, and 5 from each genome); if empty, takes full genomes
        (default empty list; list of lists of int)
    count_individuals_based_on_model -- Model object with populations and
        fitness functions used to evaluate the most fit pathogen genome in each
        host/vector in order to count only a single pathogen per host/vector, as
        opposed to all pathogens within each host/vector; if None, counts all
        pathogens (default None; None or Model)
    save_to_file -- file path and name to save model data under, no saving
        occurs if empty string (default ''; String)
    n_cores -- number of cores to parallelize processing across, if 0, all
        cores available are used (default 0; int)
    **kwargs -- additional arguents for joblib multiprocessing

    Returns:
    pandas dataframe with model sequence composition dynamics as described above
    """

    if len(populations) > 0:
        dat = cp.deepcopy( data[ data['Population'].isin( populations ) ] )
    else:
        dat = cp.deepcopy( data )

    dat = dat[ (dat['Pathogens'] != "") & (~dat['Pathogens'].isna()) ]

    if len(genomic_positions) > 0:
        print('Extracting genomic locations...')
        def extractSeq(ind):
            seqs = [ p[ genomic_positions[0]:genomic_positions[1] ]
                for p in ind['Pathogens'].split(';') ]

            return ';'.join(seqs)

        pathogens = np.array( jl.Parallel(n_jobs=n_cores, verbose=1, **kwargs) (
            jl.delayed( extractSeq ) (ind) for ind in dat
            ) )
        dat['Pathogens'] = pathogens

    if count_individuals_based_on_model is not None:
        print('Collapsing infections to individuals...')
        model = count_individuals_based_on_model
        if type_of_composition != 'Pathogens':
            raise ValueError("Computing count_individuals_based_on_model=True is only allowed for type_of_composition='Pathogens'.")

        dat['patpop'] = dat['Pathogens'].fillna('') + ':::' + dat['Population']
        unique_patpop = dat['patpop'].unique()

        for i,combination in enumerate(unique_patpop):
            print( str(i) + ' / ' + str(len(unique_patpop)) + ' combinations' )
            gen_str,pop = combination.split(':::')
            if gen_str != '':
                genomes = gen_str.split(';')
                values = np.array( [
                    model.populations[pop].fitnessHost(p)
                    for p in genomes ] )
                dominant_pathogen = genomes[ np.argmax(values) ]
                dat['patpop'] = dat['patpop'].replace(
                    combination, dominant_pathogen
                    )
            else:
                dat['patpop'] = dat['patpop'].replace( combination, '' )

        dat['Pathogens'] = dat['patpop']

    if not hosts:
        dat = dat[ dat['Organism'] != 'Host' ]
    if not vectors:
        dat = dat[ dat['Organism'] != 'Vector' ]

    all_sequences = pd.Series(
        ';'.join( dat[type_of_composition].dropna() ).split(';')
        ).str.strip()
    top_sequences = all_sequences.value_counts(ascending=False)

    if num_top_sequences < 0:
        num_top_sequences = len( top_sequences )

    if len(top_sequences) < num_top_sequences:
        num_top_sequences = len(top_sequences)

    # genomes_to_track = list(top_sequences[0:num_top_sequences].index)
    genomes_to_track = track_specific_sequences
    for genome in list(top_sequences[0:num_top_sequences].index):
        if genome not in genomes_to_track:
            genomes_to_track.append(genome)

    genome_split_data = []
    other_genome_data = pd.DataFrame(
        np.zeros( len( pd.unique( data['Time'] ) ) ),
        index=pd.unique( data['Time'] ), columns=['Other']
        )
    c = 0

    if len( ''.join(top_sequences.index) ) > 0:
        for genome in top_sequences.index:
            dat_genome = dat[ dat[type_of_composition].str.contains(
                genome, na=False, regex=False
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
                str(c) + ' / ' + str( len(top_sequences.index) )
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

    new_order = []
    for seq in track_specific_sequences:
        new_order.append(seq)
        if seq not in composition.columns:
            composition[seq] = 0

    for seq in composition.columns:
        if seq not in new_order:
            new_order.append(seq)

    composition = composition[new_order]

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
    """Create DataFrame with pairwise Hamming distances for pathogen sequences
    in data.

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
        sequences = sequences.append(
            pd.Series(track_specific_sequences)
            ).unique()

    # Fix — non parallelized
    dis_mat = np.array([[td.hamming(s1, s2) / max(len(s1),1)
        for s2 in sequences] for s1 in sequences])
            # Added the max() to avoid division by zero error triggered when
            # getPathogenDistanceHistoryDf calls this method with an empty
            # dataframe (for a timepoint with no pathogens)

    # For some reason, this code triggers a full rerun of the simulation.
    # Possible joblib bug?
    # if not n_cores:
    #     n_cores = jl.cpu_count()
    #
    # dis_mat = np.array( jl.Parallel(n_jobs=n_cores, verbose=1) (
    #     jl.delayed( lambda s1: [
    #         td.hamming(s1, s2) / max(len(s1),1) for s2 in sequences
    #         ] ) (s1) for s1 in sequences
    #     ) )

    names = sequences
    if len(seq_names) > 0:
        new_names = seq_names * int( np.ceil( len(names) / len(seq_names) ) )
        names = new_names[0:len(names)]

    dis_df = pd.DataFrame( dis_mat, index=names, columns=names )

    if len(save_to_file) > 0:
        dis_df.to_csv(save_to_file, index=True)

    return dis_df

def getPathogenDistanceHistoryDf(
        data, samples=1, num_top_sequences=-1, track_specific_sequences=[],
        seq_names=[], save_to_file="", n_cores=0):
    """Create DataFrame with pairwise Hamming distances for pathogen sequences
    in data.

    DataFrame has indexes and columns named according to genomes or argument
    seq_names, if passed. Distance is measured as percent Hamming distance from
    an optimal genome sequence.

    Arguments:
    data -- dataframe with model history as produced by saveToDf function

    Keyword arguments:
    samples -- how many timepoints to uniformly sample from the total
        timecourse; if <0, takes all timepoints (default 1; int)
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

    if samples > 0:
        samples = np.linspace(
            0, len(pd.unique(data['Time']))-1, samples
            ).astype(int)
        sampled_times = pd.unique(data['Time'])[samples]
        data = data[ data['Time'].isin(sampled_times) ]

    grouped = data.groupby('Time')
    dis_df = grouped.apply(
        lambda d: pd.melt(
            pathogenDistanceDf(
                d, num_top_sequences=num_top_sequences,
                    track_specific_sequences=track_specific_sequences,
                    seq_names=seq_names, save_to_file="", n_cores=n_cores
                    ).reset_index(),
                id_vars=['Pathogens'], var_name='Pathogen_2',
                value_name='Distance'
                )
        ).reset_index()

    dis_df.columns = ['Time','drop','Pathogen_1','Pathogen_2','Distance']
    dis_df = dis_df.drop(columns='drop')

    if len(save_to_file) > 0:
        dis_df.to_csv(save_to_file, index=False)

    return dis_df

def getGenomeTimesDf(
        data, samples=1, save_to_file="", n_cores=0, **kwargs):
    """Create DataFrame with times genomes first appeared during simulation.

    Arguments:
    data -- dataframe with model history as produced by saveToDf function

    Keyword arguments:
    samples -- how many timepoints to uniformly sample from the total
        timecourse; if <0, takes all timepoints (default 1; int)
    save_to_file -- file path and name to save model data under, no saving
        occurs if empty string (default ''; String)
    n_cores -- number of cores to parallelize across, if 0, all cores available
        are used (default 0; int)
    **kwargs -- additional arguents for joblib multiprocessing

    Returns:
    pandas dataframe with genomes and times as described above
    """

    if samples > 0:
        samples = np.linspace(
            0, len(pd.unique(data['Time']))-1, samples
            ).astype(int)
        sampled_times = pd.unique(data['Time'])[samples]
        data = data[ data['Time'].isin(sampled_times) ]

    his_dat = pd.DataFrame()

    his_dat['Sequence'] = pd.Series(
        ';'.join( dat['Pathogens'].dropna() ).split(';')
        ).str.strip()

    def getTime(seq):
        t = his_dat['Sequence'].searchsorted(seq, side='left')
        return t

    his_dat['Time_emergence'] = jl.Parallel(
        n_jobs=n_cores, verbose=1, **kwargs
        ) (
            jl.delayed( getTime ) (seq) for seq in his_dat['Sequence']
            )

    if len(save_to_file) > 0:
        dis_df.to_csv(save_to_file, index=False)

    return his_dat
