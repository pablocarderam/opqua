
"""
Simple host-host model with susceptible and infected hosts in a single
population scenario with no evolution.
"""

from opqua.model import Model # This is all you need to import to run models.

my_model = Model() # Make a new model object.
my_model.newSetup('my_setup', preset='host-host')
    # Create a new set of parameters called "my_setup" to be used to simulate
    # a population in the model. Use the default parameter set for a
    # host-host transmission model.
my_model.newPopulation('my_population', 'my_setup', num_hosts=100)
    # Create a new population of 100 hosts and 0 vectors called
    # "my_population". The population uses parameters stored in "my_setup".
my_model.addPathogensToHosts( 'my_population',{'AAAAAAAAAA':20} )
    # Add pathogens with a genome of "AAAAAAAAAA" to 20 random hosts in
    # population "my_population".
my_model.run(0,200) # Run the simulation for 200 time units.
data = my_model.saveToDataFrame('Basic_example.csv')
    # Save the model results to a table.
graph = my_model.compartmentPlot('Basic_example.png', data)
    # Plot the number of susceptible and infected hosts in the model over time.
