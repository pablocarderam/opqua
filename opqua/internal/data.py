
import numpy as np # handle arrays
import pandas as pd # data wrangling
import copy as cp
import joblib as jl

def saveToDf(history,save_to_file,n_cores=0):
    """ Saves status of model to dataframe given """

    print('Saving file...')

    if not n_cores:
        n_cores = jl.cpu_count()

    new_df = ','.join( ['Time','Population','Organism','ID','Pathogens','Protection','Alive'] ) + '\n' + \
        '\n'.join(
            jl.Parallel(n_jobs=n_cores, verbose=10) (
                jl.delayed( lambda d: ''.join(d) )
                    (
                        '\n'.join( [
                                '\n'.join( [ ','.join( [ str(time), str(pop.id), 'Host', str(host.id), '"' + ';'.join( host.pathogens.keys() ) + '"', '"' + ';'.join( host.protection_sequences ) + '"', 'True' ] )
                                    for host in pop.hosts ] ) + '\n' +
                                '\n'.join( [ ','.join( [ str(time), str(pop.id), 'Vector', str(vector.id), '"' + ';'.join( vector.pathogens.keys() ) + '"', '"' + ';'.join( vector.protection_sequences ) + '"', 'True' ] )
                                    for vector in pop.vectors ] ) + '\n' +
                                '\n'.join( [ ','.join( [ str(time), str(pop.id), 'Host', str(host.id), '"' + ';'.join( host.pathogens.keys() ) + '"', '"' + ';'.join( host.protection_sequences ) + '"', 'False' ] )
                                    for host in pop.dead_hosts ] ) + '\n' +
                                '\n'.join( [ ','.join( [ str(time), str(pop.id), 'Vector', str(vector.id), '"' + ';'.join( vector.pathogens.keys() ) + '"', '"' + ';'.join( vector.protection_sequences ) + '"', 'False' ] )
                                    for vector in pop.dead_vectors ] )
                            for id,pop in model.populations.items()
                        ] )
                    )
                for time,model in history.items()
            )
        )

    new_df = new_df.replace('\n\n','\n').replace('\n\n','\n').replace('\n\n','\n')

    file = open(save_to_file,'w')
    file.write(new_df)
    file.close()

    new_df = pd.read_csv(save_to_file)

    print('...file saved.')

    return new_df


def populationsDf(data, compartment='Infected', hosts=True, vectors=False, num_top_populations=7, track_specific_populations=[], save_to_file=""):
    ''' Populations '''

    dat = cp.deepcopy( data )

    if not hosts:
        dat = dat[ dat['Organism'] != 'Host' ]
    if not vectors:
        dat = dat[ dat['Organism'] != 'Vector' ]


    dat['Infected'] = ( dat['Pathogens'].fillna('').str.len() > 0 )
    dat['Protected'] = ( dat['Protection'].fillna('').str.len() > 0 )

    grouped = dat.groupby( [ 'Time','Population','Alive','Infected','Protected' ] ).size().reset_index(name='Number')

    compartment_names = ['Naive','Infected','Recovered','Dead']

    grouped['Compartment'] = compartment_names[3]
    grouped.loc[ ( grouped['Alive'] == True ) & ( grouped['Infected'] == False ) & ( grouped['Protected'] == False ), 'Compartment' ] = compartment_names[0]
    grouped.loc[ ( grouped['Alive'] == True ) & ( grouped['Infected'] == True ), 'Compartment' ] = compartment_names[1]
    grouped.loc[ ( grouped['Alive'] == True ) & ( grouped['Infected'] == False ) & ( grouped['Protected'] == True ), 'Compartment' ] = compartment_names[2]

    grouped = grouped[ grouped['Compartment'] == compartment ]
    grouped = grouped.drop( columns=['Alive','Infected','Protected','Compartment'] )
    grouped = grouped.pivot( columns='Population', values='Number', index='Time' ).fillna(0).reset_index('Time')

    for pop in pd.unique(data['Population']):
        if pop not in grouped.columns:
            grouped[pop] = 0


    if len(grouped.columns)-1 < num_top_populations:
        num_top_populations = len(grouped.columns)-1

    populations_to_drop = list(grouped.columns)[ num_top_populations+1: ]
    for pop in track_specific_populations:
        if pop in populations_to_drop:
            populations_to_drop = populations_to_drop.remove(pop)


    grouped['Other'] = grouped[populations_to_drop].sum(axis=1)
    if grouped['Other'].sum() == 0:
        populations_to_drop = populations_to_drop + ['Other']

    grouped = grouped.drop( columns=populations_to_drop )

    if len(save_to_file) > 0:
        grouped.to_csv(save_to_file, index=False)

    return grouped


def compartmentDf(data, populations=[], hosts=True, vectors=False, save_to_file=""):
    ''' Composition '''

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

    grouped = dat.groupby( [ 'Time','Alive','Infected','Protected' ] ).size().reset_index(name='Number')

    compartment_names = ['Naive','Infected','Recovered','Dead']

    grouped['Compartment'] = compartment_names[3]
    grouped.loc[ ( grouped['Alive'] == True ) & ( grouped['Infected'] == False ) & ( grouped['Protected'] == False ), 'Compartment' ] = compartment_names[0]
    grouped.loc[ ( grouped['Alive'] == True ) & ( grouped['Infected'] == True ), 'Compartment' ] = compartment_names[1]
    grouped.loc[ ( grouped['Alive'] == True ) & ( grouped['Infected'] == False ) & ( grouped['Protected'] == True ), 'Compartment' ] = compartment_names[2]

    grouped = grouped.drop( columns=['Alive','Infected','Protected'] )
    grouped = grouped.pivot( columns='Compartment', values='Number', index='Time' ).fillna(0).reset_index('Time')

    for comp_name in compartment_names:
        if comp_name not in grouped.columns:
            grouped[comp_name] = 0


    if len(save_to_file) > 0:
        grouped.to_csv(save_to_file, index=False)

    return grouped


def compositionDf(data, populations=[], organism='Pathogens', hosts=True, vectors=False, num_top_genomes=7, track_specific_genomes=[], save_to_file=""):
    ''' Composition '''

    if len(populations) > 0:
        dat = cp.deepcopy( data[ data['Population'].isin( populations ) ] )
    else:
        dat = cp.deepcopy( data )

    if not hosts:
        dat = dat[ dat['Organism'] != 'Host' ]
    if not vectors:
        dat = dat[ dat['Organism'] != 'Vector' ]

    all_genomes = pd.Series( ';'.join( dat[organism].dropna() ).split(';') ).str.strip()
    top_genomes = all_genomes.value_counts(ascending=False)

    if len(top_genomes) < num_top_genomes:
        num_top_genomes = len(top_genomes)

    genomes_to_track = list(top_genomes[0:num_top_genomes].index)
    for genome in track_specific_genomes:
        if genome not in genomes_to_track:
            genomes_to_track.append(genome)


    genome_split_data = []
    other_genome_data = pd.DataFrame( np.zeros( len( pd.unique( data['Time'] ) ) ), index=pd.unique( data['Time'] ), columns=['Other'] )
    c = 0

    if len( ''.join(top_genomes.index) ) > 0:
        for genome in top_genomes.index:
            dat_genome = dat[ dat[organism].str.contains(genome, na=False) ]
            grouped = dat_genome.groupby('Time').size().reset_index(name=genome)
            grouped = grouped.set_index('Time')

            if genome in genomes_to_track:
                genome_split_data.append(grouped)
            else:
                other_genome_data = other_genome_data.join( grouped, how='outer' ).fillna(0)
                other_genome_data['Other'] = other_genome_data['Other'] + other_genome_data[genome]
                other_genome_data = other_genome_data.drop(columns=[genome])

            c += 1
            print(str(c) + ' / ' + str( len(top_genomes.index) ) + ' genotypes processed.')


        if other_genome_data['Other'].sum() > 0:
            genome_split_data += [other_genome_data]
            genomes_to_track += ['Other']


    times = pd.DataFrame( np.zeros( len( pd.unique( data['Time'] ) ) ), index=pd.unique( data['Time'] ), columns=['*None*'] )
    composition = times.join( genome_split_data, how='outer' )
    composition = composition.drop(columns=['*None*']).fillna(0)
    composition = composition.reset_index()
    composition.columns = ['Time'] + list( composition.columns )[1:]

    if len(save_to_file) > 0:
        composition.to_csv(save_to_file, index=False)

    return composition


def getPathogens(dat, save_to_file=""):

    out = pd.Series( ';'.join( dat['Pathogens'].dropna() ).split(';') ).str.strip().value_counts(ascending=False).reset_index()
    out.columns = ['Pathogens'] + list( out.columns )[1:]

    if len(save_to_file) > 0:
        out.to_csv(save_to_file, index=False)

    return out

def getProtections(dat, save_to_file=""):

    out = pd.Series( ';'.join( dat['Protection'].dropna() ).split(';') ).str.strip().value_counts(ascending=False).reset_index()
    out.columns = ['Protection'] + list( out.columns )[1:]

    if len(save_to_file) > 0:
        out.to_csv(save_to_file, index=False)

    return out
