
"""
Vector-borne model with susceptible and infected hosts and vectors in a single
population scenario, showing the effect of various interventions done at
different time points.
"""

from opqua.model import Model

model = Model()
model.newSetup('my_setup',preset='vector-borne')
model.newSetup(
    'my_setup_2', contact_rate_host_vector=2e1, preset='vector-borne'
    )
    # We make a second setup with the same parameters, but duplicate the contact
    # rate.

model.newPopulation('my_population','my_setup')
model.addPathogensToHosts( 'my_population',{'AAAAAAAAAA':20} )

model.newIntervention(
    10, model.addPathogensToHosts,
    [ 'my_population', {'TTTTTTTTTT':5, 'CCCCCCCCCC':5, } ]
    )
    # At time 10, adds pathogens of genomes TTTTTTTTTT and CCCCCCCCCC to 5
    # random hosts each.

model.newIntervention( 20, model.addVectors, [ 'my_population', 10 ] )
    # At time 20, adds 10 healthy vectors to population
model.newIntervention(
    20, model.newVectorGroup,
    [ 'my_population', '10_healthy_vectors', 10, True ]
    )
    # At time 20, selects 10 healthy hosts and stores them under the group
    # ID '10_healthy_hosts'.
    # That last True argument is optional, it makes sure the random selected
    # hosts are healthy.
model.newIntervention(
    20, model.addPathogensToVectors,
    [ 'my_population', {'GGGGGGGGGG':10}, '10_healthy_vectors' ]
    )
    # At time 20, adds pathogens of genomes TTTTTTTTTT and CCCCCCCCCC to 10
    # random hosts in the '10_new_vectors' group (so, all 10 of them)
    # The last '10_healthy_vectors' argument specifies which group to sample
    # from (if not specified, sampling occurs from whole population)

model.newIntervention( 50, model.setSetup, [ 'my_population', 'my_setup_2' ] )
    # At time 50, changes the parameters of my_population to those in
    # my_setup_2.

model.newIntervention( 70, model.protectHosts, [ 'my_population', 0.75, 'A' ] )
    # At time 70, protects a random 75% of the host population from pathogens
    # with a 'A' in their genome.

output = model.run(0,100)
data = model.saveToDataFrame('Intervention_examples.csv')
graph = model.compositionPlot( 'Intervention_examples_composition.png', data )
    # Create a plot to track pathogen genotypes across time.
graph = model.compartmentPlot('Intervention_examples_compartments.png', data)
    # Also generate a normal compartment plot. Notice the total number of
    # infections in the composition plot can exceed the number of infected hosts
    # in the compartment plot. This happens because a single host infected by
    # multiple genotypes is counted twice in the former, but not the latter.
