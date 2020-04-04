
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



def compositionDf(data, populations=[], type='Pathogens', hosts=True, vectors=True, num_top_genomes=7, track_specific_genomes=[], save_to_file=""):
    ''' Composition '''

    if len(populations) > 0:
        dat = cp.deepcopy( data[ data['Population'].isin( populations ) ] )
    else:
        dat = cp.deepcopy( data )

    if not hosts:
        dat = dat[ dat['Organism'] != 'Host' ]
    if not vectors:
        dat = dat[ dat['Organism'] != 'Vector' ]

    all_genomes = pd.Series( ';'.join( dat[type].dropna() ).split(';') ).str.strip()
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
            dat_genome = dat[ dat[type].str.contains(genome, na=False) ]
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
    comp = times.join( genome_split_data, how='outer' )
    comp = comp.drop(columns=['*None*']).fillna(0)

    if len(save_to_file) > 0:
        comp.to_csv(save_to_file)

    return comp
