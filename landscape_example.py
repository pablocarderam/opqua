
"""
Here, we define a landscape with two fitness peaks, one of which is more
accessible because it has more equivalent alleles.
"""

from opqua.model import Model
import numpy as np
import scipy.stats as sp_sta
import matplotlib.pyplot as plt # plots
import seaborn as sns # pretty plots

import cProfile

my_optimal_genomes = ['AAAA','TTTT','CCCC','GGGG'] # Define an optimal genome.

# Define custom fitness function for the host:
# Fitness functions must take in 1 argument and return a positive number as a
# fitness value. Here, we take advantage of one of the preset functions, but you
# can define it any way you want!
def myHostFitness(genome):
    return np.sum([
        Model.peakLandscape(
        genome, peak_genome=my_optimal_genome, min_value=1e-1
        ) for my_optimal_genome in my_optimal_genomes
        ]) / np.sum([
            Model.peakLandscape(
            my_optimal_genomes[0], peak_genome=my_optimal_genome, min_value=1e-1
            ) for my_optimal_genome in my_optimal_genomes
            ])

model = Model()
model.newSetup( # Now, we'll define our new setup:
    'my_setup', preset='host-host', # Use default host-host parameters.
    possible_alleles='ATCG',
        # Define "letters" in the "genome", or possible alleles for each locus.
        # Each locus can have different possible alleles if you define this
        # argument as a list of strings, but here, we take the simplest
        # approach.
    allele_groups_host=[['A','TCG']],
        # available alleles are grouped in two groups of alleles with equivalent
        # behavior: B, C, and D are phenotypically identical across all loci
    num_loci=len(my_optimal_genomes[0]),
        # Define length of "genome", or total number of alleles.
    fitnessHost=myHostFitness,
        # Assign the fitness function we created (could be a lambda function)
    mutate_in_host=1e-2
        # Modify de novo mutation rate of pathogens when in host to get some
        # evolution!
    )

model.newLandscape('my_setup', 'my_landscape', fitnessFunc=myHostFitness,
          population_threshold=1e10, max_depth=4)
    # defines new fitness landscape object based on setup parameters; depth=4
    # and a large population threshold means it will be fully evaluated

gens = 2000
s = 0.01
land = model.setups['my_setup'].landscapes['my_landscape']
sur_mut = land.survivalProbabilities(
    s, gens, land.population_size
    )
sur_neu = land.survivalProbabilities(
    0, gens, land.population_size
    )
mut_wait_time = np.array([
    sp_sta.geom.pmf( x+1, land.mutate )
    for x in range( gens )
    ])
mult_mut = np.multiply(sur_mut,mut_wait_time).cumsum()
mult_neu = np.multiply(sur_neu,mut_wait_time).cumsum()
quo = np.divide( mult_mut, mult_neu )
# print(sur)
# print(mut_wait_time)
# print(mult)
print(sur_mut.sum(),sur_neu.sum(),mut_wait_time.sum(),mult_mut[-1],mult_neu[-1],quo[-1])
ax = plt.subplot(6, 1, 1)
ax.plot(sur_mut)
ax = plt.subplot(6, 1, 2)
ax.plot(sur_neu)
ax = plt.subplot(6, 1, 3)
ax.plot(mut_wait_time)
ax = plt.subplot(6, 1, 4)
ax.plot(mult_mut)
ax = plt.subplot(6, 1, 5)
ax.plot(mult_neu)
ax = plt.subplot(6, 1, 6)
ax.plot(quo)
plt.savefig('surv.png', bbox_inches='tight')

cProfile.run("model.mapLandscape('my_setup', 'my_landscape', 'AAAA')",sort='cumtime')
# model.mapLandscape('my_setup', 'my_landscape', 'AAAA')
    # evaluates and maps mutation establishment rates across landscape

model.saveLandscape('my_setup', 'my_landscape', 'landscape.csv')
model.loadLandscape('my_setup', 'my_landscape', 'landscape.csv')
print(model.setups['my_setup'].landscapes['my_landscape'].mutation_network)

model.visualizeMutationNetwork('my_setup', 'my_landscape', 'mutation_network')
    # creates interactive visualization
