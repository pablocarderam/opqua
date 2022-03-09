
"""
Host-host transmission model with susceptible and infected hosts in a single
population scenario, illustrating pathogen evolution through de novo mutations
and acquired immunity.
"""

from opqua.model import Model

# Define custom immunity function for the host:
# Immunity functions must take in 2 arguments and return a 0-1 immunity
# coefficient value. Here, we take advantage of one of the preset functions,
# but you can define it any way you want!
def myImmuneWeights(genome,immune_seq):
    return Model.perfectMatchImmunity(
        genome, immune_seq, weight=1
        )
        # Perfect matching: any mismatches between pathogen genome and immune
        # memory generate no immunity. Only perfect matches generate 100%
        # immunity.

model = Model()
model.newSetup( # Now, we'll define our new setup:
    'my_setup', preset='host-host', # Use default host-host parameters.
    possible_alleles='AB',
        # Define "letters" in the "genome", or possible alleles for each locus.
        # Each locus can have different possible alleles if you define this
        # argument as a list of strings, but here, we take the simplest
        # approach.
    num_loci=3,
        # Define length of "genome", or total number of alleles.
    mutate_in_host=1e-2,
        # Modify de novo mutation rate of pathogens when in host to get some
        # evolution!
    immunity_acquisition_rate_host=1e-2,
        # rate at which immunity is acquired within infected individuals
    immunity_loss_rate_host=0,
        # rate at which immunity is lost within infected individuals
    immunityWeightsHost=myImmuneWeights,
        # immunity function that evaluates effect of immune memory on pathogen
        # genomes.
    )

model.newPopulation('my_population','my_setup', num_hosts=100)
model.addPathogensToHosts( 'my_population',{'AAA':10} )
    # We will start off the simulation with a suboptimal pathogen genome,
    # "AAA". Throughout the course of the simulation, we should see this genome
    # be outcompeted by other pathogen genotypes as the host population acquires
    # resistance to each genome.
model.run(0,1000,time_sampling=0)
data = model.saveToDataFrame(
    'immunity_example.csv'
    )

graph_composition = model.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'immunity_example_pathogen_composition.png', data,
    type_of_composition='Pathogens',
    num_top_sequences=8,
        # Track the 8 most represented genomes overall (only 8 possible).
    )

graph_composition = model.compositionPlot(
        # Create a plot of genotypes in the hosts' immune memories across time.
    'immunity_example_immunity_composition.png', data,
    type_of_composition='Immunity', y_label='Genomes in immune memory',
    num_top_sequences=8,
        # Track the 8 most represented genomes overall (only 8 possible).
    )

graph_compartments = model.compartmentPlot(
    'immunity_example_reassortment_compartments.png', data
    )
    # Also generate a normal compartment plot. Notice the total number of
    # infections in the composition plot can exceed the number of infected hosts
    # in the compartment plot. This happens because a single host infected by
    # multiple genotypes is counted twice in the former, but not the latter.
