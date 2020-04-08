
'''
Example with a model
'''

from opqua.model import Model

my_optimal_genome = 'BEST' # define an optimal genome

# Define custom fitness function for the host:
# Fitness functions must take in 1 argument and return a positive number
# as a fitness value. Here, we take advantage of one of the preset
# functions, but you can define it any way you want!
def myHostFitness(genome):
    return Model.stabilizingSelection(
        genome, optimal_genome=my_optimal_genome, min_fitness=1e-10
        )
        # Stabilizing selection: any deviation from the "optimal genome"
        # sequence results in an exponential decay in fitness to the min_fitness
        # value at the maximum possible distance.

model = Model()
model.newSetup(
    'my_setup', preset='host-host', # use default host-host parameters
    possible_alleles='ABDEST', # letters in the "genome"
    num_loci=len(my_optimal_genome), # length of "genome"
    fitnessHost=myHostFitness, # feed in custom
    mutate_in_host=5e-2
    )
    # uses default parameters except for possible alleles per locus and host
    # fitness function

model.newPopulation('my_population','my_setup')
model.addPathogensToHosts( 'my_population',{'BADD':10} )
model.run(0,200)
data = model.saveToDataFrame('Stabilizing_selection.csv')
graph_composition = model.compositionPlot(
    'Stabilizing_selection_composition.png', data, num_top_sequences=6,
    track_specific_genomes=['BADD']
    )
graph_compartments = model.compartmentPlot(
    'Stabilizing_selection_compartments.png', data
    )
