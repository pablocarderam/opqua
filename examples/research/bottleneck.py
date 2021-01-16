
"""
Host-host transmission model with susceptible and infected hosts in a single
population scenario, illustrating pathogen evolution through a simple
bottleneck.
"""

from opqua.model import Model

model = Model()
model.newSetup( # Now, we'll define our new setup:
    'my_setup', preset='host-host', # Use default host-host parameters.
    num_loci=100,
        # Define length of "genome", or total number of alleles.
    # recovery_rate_host = 2e-1,
    # contact_rate_host_host = 1e2,
    mutate_in_host=1e-1 # Modify mutation rate to get some evolution!
    )

model.newPopulation('my_population','my_setup')
model.addPathogensToHosts( 'my_population',{'T'*100:10} )
    # We will start off the simulation with 10 cases of a 100-mer 'T'
    # homopolymer genome

model.newIntervention( 100, model.protectHosts,
                      [ 'my_population', 0.8, 'T'*30 ] )
    # At time 100, vaccinates random 80% of hosts against genomes with 10-mer Ts

model.run(0,200)
data = model.saveToDataFrame('Bottleneck.csv')

graph_composition = model.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'Bottleneck_composition.png', data,
    num_top_sequences=6,
        # Track the 6 most represented genomes overall (remaining genotypes are
        # lumped into the "Other" category).
    )

model.pathogenDistanceHistory(data,samples=100,
    save_to_file='pairwise_distance_history.csv')

graph_clustermap = model.clustermap(
    'Bottleneck_clustermap.png', data,
    save_data_to_file='Bottleneck_pairwise_distances.csv',
    num_top_sequences=15
    )
    # Generate a heatmap and dendrogram for the top 15 genomes.
    # Besides creating the plot,
    # outputs the pairwise distance matrix to a csv file as well.

graph_compartments = model.compartmentPlot(
    'Bottleneck_compartments.png', data
    )
    # Also generate a normal compartment plot. Notice the total number of
    # infections in the composition plot can exceed the number of infected hosts
    # in the compartment plot. This happens because a single host infected by
    # multiple genotypes is counted twice in the former, but not the latter.
