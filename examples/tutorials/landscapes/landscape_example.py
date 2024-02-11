
"""
Here, we define a landscape with two fitness peaks, one of which is more
accessible because it has more equivalent alleles.
"""

from opqua.model import Model
import numpy as np

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

model.mapLandscape('my_setup', 'my_landscape', 'AAAA')
    # evaluates and maps mutation establishment rates across landscape

model.visualizeMutationNetwork('my_setup', 'my_landscape', 'mutation_network')
    # creates interactive visualization
