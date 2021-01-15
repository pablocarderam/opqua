"""One single population in which multiple genotypes are coexisting sympatrically
-Low transmission (1)
-Fitnes function sets LLLLL as the optimal genotype and the other genotypes' fitness is given by the hammond distance
to the optimal genome.
-Number of loci = 20
-Possible alleles = 4
-Mutation in host is set high (1e-1)
Note: If the contact rate is <=1, genotypes dissapear after aproximmately 50 time units"""

from opqua.model import Model
import numpy as np
import textdistance as td

my_model = Model()

def fitnessGenotype(genome):
    WT = 'A' * len(genome)
    Resistant = 'M' * len(genome)
    if genome == WT or genome==Resistant:
        return 1
    distanceWT = td.hamming(Resistant, genome) /len(genome)
    distanceRes = td.hamming(WT, genome) / len(genome)
    distance = max(distanceWT, distanceRes)
    fitness = np.exp(np.log(1e-10) * distance)
    return fitness

my_model.newSetup('my_setup',preset='vector-borne',num_loci=20, possible_alleles='AMDG',contact_rate_host_vector=5,
                  fitnessHost=fitnessGenotype, mutate_in_host=1e-4, recombine_in_vector=1)

my_model.newPopulation('my_population', 'my_setup', num_hosts=500, num_vectors=500)
    # Create a new population of 100 hosts and 100 vectors called
    # "my_population". The population uses parameters stored in "my_setup".
my_model.addPathogensToHosts('my_population',{('D'*20):50,('M'*20):50})
my_model.run(0,200) # Run the simulation for 200 time units.
data = my_model.saveToDataFrame('LowFitnessLGD.csv')
graph = my_model.compositionPlot("LowFitLowGenotypes.png", data)
graph1 = my_model.compartmentPlot('LowFit_Hosts.png', data)
