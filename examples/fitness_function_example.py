
"""
Host-host transmission model with susceptible and infected hosts in a single
population scenario, illustrating pathogen evolution through intra-host
competition.

When two pathogens with different genomes meet in the same host (or vector),
the pathogen with the most fit genome has a higher probability of being
transmitted to anothre host (or vector). Here, we define a landscape of
stabilizing selection where there is an optimal genome and every other genome is
less fit, but fitness functions can be defined in any arbitrary way (accounting
for multiple peaks, for instance, or special cases of a genome).
"""

from opqua.model import Model

my_optimal_genome = 'BEST' # Define an optimal genome.

# Define custom fitness function for the host:
# Fitness functions must take in 1 argument and return a positive number as a
# fitness value. Here, we take advantage of one of the preset functions, but you
# can define it any way you want!
def myHostFitness(genome):
    return Model.stabilizingSelection(
        genome, optimal_genome=my_optimal_genome, min_fitness=1e-10
        )
        # Stabilizing selection: any deviation from the "optimal genome"
        # sequence results in an exponential decay in fitness to the min_fitness
        # value at the maximum possible distance. Here we use strong selection,
        # with a very low minimum fitness.

model = Model()
model.newSetup( # Now, we'll define our new setup:
    'my_setup', preset='host-host', # Use default host-host parameters.
    possible_alleles='ABDEST',
        # Define "letters" in the "genome", or possible alleles for each locus.
        # Each locus can have different possible alleles if you define this
        # argument as a list of strings, but here, we take the simplest
        # approach.
    num_loci=len(my_optimal_genome),
        # Define length of "genome", or total number of alleles.
    fitnessHost=myHostFitness,
        # Assign the fitness function we created (could be a lambda function)
    mutate_in_host=5e-2 # Modify mutation rate to get some evolution!
    )

model.newPopulation('my_population','my_setup')
model.addPathogensToHosts( 'my_population',{'BADD':10} )
    # We will start off the simulation with a suboptimal pathogen genome,
    # "BADD". Throughout the course of the simulation, we should see this genome
    # be outcompeted by more optimal pathogen genotypes, culminating in the
    # optimal genome, "BEST", which outcompetes all others.
model.run(0,200)
data = model.saveToDataFrame('Stabilizing_selection.csv')

graph_composition = model.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'Stabilizing_selection_composition.png', data,
    num_top_sequences=6,
        # Track the 6 most represented genomes overall (remaining genotypes are
        # lumped into the "Other" category).
    track_specific_genomes=['BADD']
        # Include the initial genome in the graph if it isn't in the top 6.
    )

graph_compartments = model.compartmentPlot(
    'Stabilizing_selection_compartments.png', data
    )
    # Also generate a normal compartment plot. Notice the total number of
    # infections in the composition plot can exceed the number of infected hosts
    # in the compartment plot. This happens because a single host infected by
    # multiple genotypes is counted twice in the former, but not the latter.
