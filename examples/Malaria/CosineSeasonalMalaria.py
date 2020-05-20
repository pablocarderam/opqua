"""Seasonal malaria simulation with single population of infected/susceptible hosts and vectors
-Standard parameters (default setup)
-Two whole periods of the absolute value of cosine function is given as the transmission rate"""

from opqua.model import Model
import numpy as np
import textdistance as td
my_model = Model() # Make a new model object.

def fitnessGenotype(genome):
    WT = 'A' * len(genome)
    Resistant = 'Z' * len(genome)
    if genome == WT or genome==Resistant:
        return 1
    distanceWT = td.hamming(Resistant, genome) /len(genome)
    distanceRes = td.hamming(WT, genome) / len(genome)
    distance = max(distanceWT, distanceRes)
    fitness = np.exp(np.log(1e-10) * distance)
    return fitness

my_model.newSetup('my_setup',preset='vector-borne', fitnessHost=fitnessGenotype, possible_alleles="AFMNUZ")
    # Create a new set of parameters called "my_setup" to be used to simulate
    # a population in the model. Use the default parameter set for a
    # vector-borne model.
my_model.newPopulation('my_population', 'my_setup', num_hosts=100, num_vectors=100)
    # Create a new population of 100 hosts and 100 vectors called
    # "my_population". The population uses parameters stored in "my_setup".
my_model.addPathogensToHosts( 'my_population',{('M'*10):5,('N'*10):5,('U'*10):5,('F'*10):5})
    # Add pathogens with a genome of "AAAAAAAAAA" to 20 random hosts in
    # population "my_population".
for i in range(0, 190, 5):
    rate = 100 * np.abs((np.cos(np.pi * i / 50))) + 1
    my_model.newSetup("i_setup", preset='vector-borne', contact_rate_host_vector=rate)
    my_model.newIntervention(i, my_model.setSetup, ['my_population', 'i_setup'])

my_model.run(0,200) # Run the simulation for 200 time units.
data = my_model.saveToDataFrame('Cosine_Seasonal.csv')
    # Save the model results to a table.
graph = my_model.compartmentPlot('Cosine_Seasonal_Hosts.png', data)
#Plot the infected and susceptible hosts across time
graph2 = my_model.compositionPlot("Cosine_Seasonal_Genotypes.png",data)
#Plot the different genotypes across time