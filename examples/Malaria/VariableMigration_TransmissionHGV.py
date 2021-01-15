"""Three different populations with a variable rate of migration between them
-Each population has a different transmission rate (contact_rate_host_vector)
-Population1: 1e3 , Population2: 1e1, Population3: 1e-1
-High genetic variation: 5 possible alleles
-Mutation in vector is given a value of 1e-4 for all populations
-Migration among populations increases with time
-The files for this simulation are saved with HGV at the end, which stands for high genetic variation. This is done to set
apart this simulation from another with same parameters but low genetic variation"""

import numpy as np
from opqua.model import Model # This is all you need to import to run models.

my_model = Model() # Make a new model object.

my_model.newSetup('highT',preset='vector-borne', num_loci=5,possible_alleles="DRSTW",contact_rate_host_vector=1e3, mutate_in_host=1e-4)
my_model.newSetup('mediumT',preset='vector-borne', num_loci=5,possible_alleles="DRSTW",contact_rate_host_vector=1e1,mutate_in_host=1e-4)
my_model.newSetup('lowT',preset='vector-borne', num_loci=5,possible_alleles="DRSTW",contact_rate_host_vector=1e-1, mutate_in_host=1e-4)
    # Create a new set of parameters called "my_setup" to be used to simulate
    # a population in the model. Use the default parameter set for a
    # vector-borne model.

my_model.newPopulation('population_highT', 'highT', num_hosts=30, num_vectors=30)
my_model.newPopulation('population_mediumT', 'mediumT', num_hosts=30, num_vectors=30)
my_model.newPopulation('population_lowT', 'lowT', num_hosts=30, num_vectors=30)

my_model.addPathogensToHosts( 'population_highT',{'DDDDD':10})
my_model.addPathogensToHosts( 'population_mediumT',{'SSSSS':10})
my_model.addPathogensToHosts( 'population_lowT',{'WWWWW':10})

for i in range(0,200,1):
    migration=0.01*np.sqrt(i)
    my_model.linkPopulations('population_highT','population_mediumT',migration)
    my_model.linkPopulations('population_mediumT','population_highT',migration)
    my_model.linkPopulations('population_highT','population_lowT',migration)
    my_model.linkPopulations('population_lowT','population_highT',migration)
    my_model.linkPopulations('population_lowT','population_mediumT',migration)
    my_model.linkPopulations('population_mediumT','population_lowT',migration)
my_model.run(0,200) # Run the simulation for 200 time units.
data = my_model.saveToDataFrame('VarMigr_TransmHGV.csv')
    # Save the model results to a table.
graph = my_model.compartmentPlot('VarMigr_Transm_HostsHGV.png', data)
    # Plot the number of susceptible and infected hosts in the model over time.
graph2 = my_model.compositionPlot('VarMigr_Transm_GenotypesHGV.png',data)
    #Plot the track genotypes across time
graph3= my_model.populationsPlot('VarMigr_Transm_PopulationHGV',data)
    #Track the number of infected per population over time
