"""One single population in which multiple genotypes are coexisting sympatrically.
-Low transmission (1)
-Fitnes function sets wild type and resistant genotypes as the optimal genotypes and the other genotypes' fitness is
given by the max hammond distance to either of the optimal genotypes.
-Number of loci = 20
-Possible alleles = 8"""

from opqua.model import Model
import numpy as np
import textdistance as td

my_model = Model()

def myHostFitness(genome):
    WT = 'A' * len(genome)
    Resistant = 'M' * len(genome)
    if genome == WT or genome == Resistant:
        return 1
    distanceWT = td.hamming(Resistant, genome) / len(genome)
    distanceRes = td.hamming(WT, genome) / len(genome)
    distance = max(distanceWT, distanceRes)
    fitness = np.exp(np.log(1e-10) * distance)
    return fitness
my_model.newSetup('my_setup',preset='vector-borne',num_loci=20, possible_alleles='ABDGEJLM',contact_rate_host_vector=5,
                  fitnessHost=myHostFitness, mutate_in_host=1e-4, recombine_in_vector=1)

my_model.newPopulation('my_population', 'my_setup', num_hosts=500, num_vectors=500)
    # Create a new population of 500 hosts and 500 vectors called
    # "my_population". The population uses parameters stored in "my_setup".
my_model.addPathogensToHosts('my_population',{('B'*10)+('L'*10):30,('D'*10)+('J'*10):30,('G'*10)+('E'*10):30})
my_model.run(0,200) # Run the simulation for 200 time units.
data = my_model.saveToDataFrame('LowFitnessHGD.csv')
graph = my_model.compositionPlot("LowFitHighGenotypes.png", data)
graph1 = my_model.compartmentPlot('LowFit_HostsHGD.png', data)