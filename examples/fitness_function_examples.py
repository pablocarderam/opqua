
'''
Example with a model
'''

from opqua.model import Model

# Used for custom fitness function:
import numpy as np

# Define custom fitness function for the host
my_optimal_genome = 'BEST'

def myHostFitness(genome):
    """ Fitness functions must take in 1 argument and retuen a positive number
        as a fitness value. Here, we take advantage of one of the presets,
        but you can define it any way you want! """

    return Model.stabilizingSelection( genome, optimal_genome=my_optimal_genome, min_fitness=1e-5)
        # Stabilizing selection: any deviation from the "optimal genome" sequence
        # results in an exponential decay in fitness to the min_fitness value at
        # the maximum possible distance.

my_model = Model()
my_model.newSetup( 'my_setup', default='host-host', # use default host-host parameters
                    possible_alleles='ABDEST', # letters in the "genome"
                    num_loci=len(my_optimal_genome), # length of "genome"
                    fitnessHost=myHostFitness, # feed in custom
                    mutate_in_host=0.1
                    ) # uses default parameters except for possible alleles per locus and host fitness function

my_model.newPopulation('my_population','my_setup')
my_model.addPathogensToHosts( 'my_population',{'BADD':4} )
output = my_model.run(0,100,'StabilizingSelection.csv')
