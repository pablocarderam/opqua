"""One population of 100 hosts and 100 vectors in which the fitness of the genotype is given by the maximum hamming
distance to the wild type genotype or the resistant genotype (which have the greatest fitness values)
-The number of loci is 5, with 5 possible alleles
-The fitness function is given by the maximum hamming distance to either the wild type or the resistant genotypes,
which have the greatest fitness (1)
-Default setup"""
from opqua.model import Model
import textdistance as td
import numpy as np

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

my_model = Model() # Make a new model object.
my_model.newSetup('my_setup',preset='vector-borne', num_loci=5,possible_alleles="AMNPZ",
                  fitnessHost=fitnessGenotype, recombine_in_vector=1)
    # Create a new set of parameters called "my_setup" to be used to simulate
    # a population in the model. Use the default parameter set for a
    # vector-borne model.
my_model.newPopulation(
    'my_population', 'my_setup', num_hosts=100, num_vectors=100)
    # Create a new population of 100 hosts and 100 vectors called
    # "my_population". The population uses parameters stored in "my_setup".
my_model.addPathogensToHosts( 'my_population',{('M'*5):10,('P'*5):10,('N'*5):10} )
    # Add pathogens with a genome of "MMMMM","PPPPP", "NNNNN", to 10 random hosts for each genotype in
    # population "my_population".
my_model.run(0,200) # Run the simulation for 200 time units.
data = my_model.saveToDataFrame('FitnessValley.csv')
    # Save the model results to a table.
graph_2= my_model.compositionPlot('FitnessValleyComp.png', data, num_top_sequences=10,
                                  track_specific_sequences=['ZZZZZ','AAAAA'])
#Plot the different genotypes across time, tracking the wild type
# and the resistant ones (AAAAA and ZZZZZ respectively)