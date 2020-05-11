"""Two isolated populations with a fitness function in which the optimal genomes are the wild type and the resistant.
This fitness function is used for both fitnesHost and fitnessVector.
The current
For both populations the number of hosts and vectors is 100, and the default setup is used
For the first population (later referenced as popLT) the following parameters are changed:
-Contact_rate_host_vector = 5e-1
For the second population (later referenced as popHT) the following parameters are given:
-Contact_rate_host_vector = 5e1
Note:
-If recombine or mutation are not sufficiently high, there is little to no apearance of
different genotypes than the originally inserted with the addPathogensToHost function.
This results in the dissapearance of infected hosts at small times when contact_rat_hosts_vectors
is sufficiently small (usually less than 1e-1)."""

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

my_model.newSetup('my_setup',preset='vector-borne',num_loci=5, possible_alleles='AMDG',contact_rate_host_vector=5e-1,
                 fitnessHost=fitnessGenotype, fitnessVector=fitnessGenotype)
my_model.newSetup('my_setup2',preset='vector-borne',num_loci=5, possible_alleles='AMDG',contact_rate_host_vector=5e1,
                 fitnessHost=fitnessGenotype, fitnessVector=fitnessGenotype)
my_model.newPopulation('popLT', 'my_setup', num_hosts=100, num_vectors=100)
my_model.newPopulation('popHT', 'my_setup2', num_hosts=100, num_vectors=100)
# Create two new populations of 100 hosts and 100 vectors called
# "my_population" and "my_population2. The populations use parameters stored in "my_setup" and "my_setup2".

my_model.addPathogensToHosts('popLT',{('D'*5):10,('G'*5):10,'MMMMM':10})
my_model.addPathogensToHosts('popHT',{('D'*5):10,('G'*5):10,'MMMMM':10})
#To both populations are added 20 pathogens with non optimal genotypes

for i in range(5,200,3):
    my_model.newIntervention( i , my_model.treatHosts ,
                             [ "popLT" , 0.1 , ['M'*5]]
                             )
    my_model.newIntervention(i, my_model.treatHosts,
                             ["popHT", 0.1, ['M' * 5]]
                             )
    #Treat 10% of the infected hosts every 5 time units in both populations

my_model.run(0,200) # Run the simulation for 200 time units.
data = my_model.saveToDataFrame('DrugIntervention.csv')

graph11 = my_model.compositionPlot("DrugInterventionGenotypesLT.png", data, populations=["popLT"],
                                   track_specific_sequences=['M'*5])
graph12 = my_model.compositionPlot("DrugInterventionGenotypesHT.png", data, populations=["popHT"],
                                   track_specific_sequences=['M'*5])
graph2 = my_model.populationsPlot("DrugIntervention_Populations.png",data)