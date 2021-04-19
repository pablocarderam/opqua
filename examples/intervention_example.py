
"""
Vector-borne model with susceptible and infected hosts and vectors in a single
population scenario, showing the effect of various interventions done at
different time points.
"""

from opqua.model import Model

model = Model()
model.newSetup('my_setup',preset='vector-borne')
model.newSetup(
    'my_setup_2', contact_rate_host_vector=4e-1, preset='vector-borne'
    )
    # We make a second setup with the same parameters, but duplicate the contact
    # rate.

model.newPopulation('my_population','my_setup')
model.addPathogensToHosts( 'my_population',{'AAAAAAAAAA':20} )

model.newIntervention(
    20, model.addPathogensToHosts,
    [ 'my_population', {'TTTTTTTTTT':5, 'CCCCCCCCCC':5, } ]
    )
    # At time 20, adds pathogens of genomes TTTTTTTTTT and CCCCCCCCCC to 5
    # random hosts each.

model.newIntervention( 50, model.addVectors, [ 'my_population', 10 ] )
    # At time 50, adds 10 healthy vectors to population
model.newIntervention(
    50, model.newVectorGroup,
    [ 'my_population', '10_healthy_vectors', 10, 'healthy' ]
    )
    # At time 50, selects 10 healthy hosts and stores them under the group
    # ID '10_healthy_hosts'.
model.newIntervention(
    50, model.addPathogensToVectors,
    [ 'my_population', {'GGGGGGGGGG':10}, '10_healthy_vectors' ]
    )
    # At time 50, adds pathogens of genomes GGGGGGGGGG to 10
    # random hosts in the '10_new_vectors' group (so, all 10 of them)
    # The last '10_healthy_vectors' argument specifies which group to sample
    # from (if not specified, sampling occurs from whole population)

model.newIntervention( 100, model.setSetup, [ 'my_population', 'my_setup_2' ] )
    # At time 100, changes the parameters of my_population to those in
    # my_setup_2, with twice the contact rate

model.newIntervention(
    150, model.newHostGroup,
    [ 'my_population', 'treated_hosts', -1, 'infected' ]
    )
    # At time 150, selects 100% of infected hosts and stores them under the
    # group ID 'treated_hosts'.
model.newIntervention(
    150, model.newVectorGroup,
    [ 'my_population', 'treated_vectors', -1, 'infected' ]
    )
    # At time 150, selects 100% of infected vectors and stores them under the
    # group ID 'treated_vectors'.
model.newIntervention( 150, model.treatHosts, [ 'my_population', 1, 'GGGGGGGGGG', 'treated_hosts' ] )
    # At time 50, treat 100% of the "treated_hosts" population with a treatment
    # that kills pathogens unless they contain a 'GGGGGGGGGG' sequence in their
    # genome.
model.newIntervention( 150, model.treatVectors, [ 'my_population', 1, 'GGGGGGGGGG', 'treated_vectors' ] )
    # At time 50, treat 100% of the "treated_vectors" population with a
    # treatment that kills pathogens unless they contain a 'GGGGGGGGGG' sequence
    # in their genome.

model.newIntervention(
    250, model.newHostGroup,
    [ 'my_population', 'vaccinated', 0.85, 'any' ]
    )
    # At time 250, selects 85% of random hosts and stores them under the group
    # ID 'vaccinated'. They may be healthy or infected.
model.newIntervention( 250, model.protectHosts, [ 'my_population', 1, 'GGGGGGGGGG', 'vaccinated' ] )
# model.newIntervention( 250, model.protectHosts, [ 'my_population', 1, 'C', 'vaccinated' ] )
# model.newIntervention( 250, model.protectHosts, [ 'my_population', 1, 'T', 'vaccinated' ] )
    # At time 250, protects 100% of the vaccinated group from pathogens
    # with a 'GGGGGGGGGG' sequence in their genome.

output = model.run(0,400)
data = model.saveToDataFrame('Intervention_examples.csv')
graph = model.compositionPlot( 'Intervention_examples_composition.png', data )
    # Create a plot to track pathogen genotypes across time.
graph = model.compartmentPlot('Intervention_examples_compartments.png', data)
    # Also generate a normal compartment plot. Notice the total number of
    # infections in the composition plot can exceed the number of infected hosts
    # in the compartment plot. This happens because a single host infected by
    # multiple genotypes is counted twice in the former, but not the latter.
