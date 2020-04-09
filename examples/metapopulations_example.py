
"""
Vector-borne model with susceptible and infected hosts and vectors, showing a
metapopulation model setup with multiple populations connected to each other.
"""

from opqua.model import Model

model = Model()
model.newSetup('setup_normal', preset='vector-borne')
model.newSetup(
    'setup_cluster',
    contact_rate_host_vector = ( 2 *
        model.setups['setup_normal'].contact_rate_host_vector ),
    preset='vector-borne'
    ) # uses default parameters but doubles contact rate of the first setup

model.newPopulation('population_A','setup_normal', num_hosts=20, num_vectors=20)
model.newPopulation('population_B','setup_normal', num_hosts=20, num_vectors=20)
    # Create two populations that will be connected.
model.newPopulation(
    'isolated_population','setup_normal', num_hosts=20, num_vectors=20
    ) # A third population will remain isolated.

model.createInterconnectedPopulations(
    5,1e-2,'clustered_population_','setup_cluster', num_hosts=20, num_vectors=20
    )
    # Create a cluster of 5 populations connected to each other with a migration
    # rate of 1e-2 between each of them in both directions. Each population has
    # an numbered ID with the prefix "clustered_population_", has the parameters
    # defined in the "setup_cluster" setup, and has 20 hosts and vectors.
model.linkPopulations('population_A','clustered_population_4',1e-2)
    # We link population_A to one of the clustered populations with a one-way
    # migration rate of 1e-2.
model.linkPopulations('population_A','population_B',1e-2)
    # We link population_A to population_B with a one-way migration rate of
    # 1e-2.

model.addPathogensToHosts( 'population_A',{'AAAAAAAAAA':5} )
    # population_A starts with AAAAAAAAAA genotype pathogens.
model.addPathogensToHosts( 'population_B',{'GGGGGGGGGG':5} )
    # population_B starts with GGGGGGGGGG genotype pathogens.

output = model.run(0,100)
data = model.saveToDataFrame('Metapopulations_example.csv')
graph = model.populationsPlot( # Plot infected hosts per population over time.
    'Metapopulations_example.png', data,
    track_specific_populations=['isolated_population'],
        # Make sure to plot th isolated population totals if not in the top
        # infected populations.
    x_label='Infected hosts' # change x label
    )
