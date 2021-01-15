"""Three different populations, tow of which are interconnected at an increasing rate of migration.
 The third population is isolated but increases its standing variation across time
-The two firsts populations have a different genetic variation, the third one has a genetic variation that increases
over time.
-Population1: 5 possible alleles, Population2: 3 possible alleles, Population3: initially 1 possible allele. Every 30
time units, 3 possible alleles are added and 2 pathogens with each of those genotypes are added to the vectors.
-Low transmission (contact_rate_host_vector) is constant for all populations: 1
-Mutation in vector is given a value of 1 for all populations
-Number of loci is set to 5
-The fittness for a genotype, both in host and vector, is given by a function that sets the wild type and the resistant
genotypes as the optimum genotypes. The fitness of other genotypes is then calculated by taking the maximum hamming
distance to the optimal genotypes."""

from opqua.model import Model # This is all you need to import to run models.
import numpy as np
import textdistance as td

def fitnessGenotype(genome):
    WT = 'A' * len(genome)
    Resistant = 'M' * len(genome)
    if genome == WT or genome==Resistant:
        return 1
    distanceWT = td.hamming(Resistant, genome) /len(genome)
    distanceRes = td.hamming(WT, genome) / len(genome)
    distance = max(distanceWT, distanceRes)
    fitness = np.exp(np.log(1e-5) * distance)
    return fitness

my_model = Model() # Make a new model object.

my_model.newSetup('highGV_LT',preset='vector-borne', num_loci=5,possible_alleles="DRSTW",
                  contact_rate_host_vector=5, mutate_in_host=1e-2, fitnessVector=fitnessGenotype,
                  fitnessHost=fitnessGenotype)
my_model.newSetup('mediumGV_LT',preset='vector-borne', num_loci=5,possible_alleles="BCY",
                  contact_rate_host_vector=5,mutate_in_host=1e-2, fitnessVector=fitnessGenotype,
                  fitnessHost=fitnessGenotype)
my_model.newSetup('isolatedGV_LT',preset='vector-borne', num_loci=5,possible_alleles="A",
                  contact_rate_host_vector=1, mutate_in_host=1e-2, fitnessVector=fitnessGenotype,
                  fitnessHost=fitnessGenotype)
    # Create three new set of parameters to be used to simulate
    # the three populations in the model. Use the default parameter set for a
    # vector-borne model.

my_model.newPopulation('population_highGV', 'highGV_LT', num_hosts=50, num_vectors=50)
my_model.newPopulation('population_mediumGV', 'mediumGV_LT', num_hosts=50, num_vectors=50)
my_model.newPopulation('isolatedGV', 'isolatedGV_LT', num_hosts=50, num_vectors=50)

my_model.addPathogensToHosts( 'population_highGV',{'DRSSW':2})
my_model.addPathogensToHosts( 'population_mediumGV',{'BBCYY':2})
my_model.addPathogensToHosts( "isolatedGV",{'AAAAA':2})

for i in range(0,200,1):
    migration=0.01*np.sqrt(i)
    my_model.linkPopulations('population_highGV','population_mediumGV',migration)
    my_model.linkPopulations('population_mediumGV','population_highGV',migration)
#Populations 1 and 2 (with medium and high genetic variation, are linked with an increasing transmission rate

for time in range(30, 200, 30):
    letters = 'ABCDEFGHIJKLM' #String of all possible genotypes
    i = int(time / 30) #Convert time into indexes that go from 1 to 6
    pos_al = letters[:(2 * i) + 1] #Use the index obtained above to slice the "letters" string. This way every
    # intervention you can add two more possible alleles

    my_model.newSetup('i_setup', preset='vector-borne', possible_alleles=pos_al, num_loci=5, contact_rate_host_vector=5,
                      mutate_in_host=1e-2, fitnessVector=fitnessGenotype, fitnessHost=fitnessGenotype)
    my_model.newIntervention(time, my_model.setSetup, ['isolatedGV', 'i_setup'])
    my_model.newIntervention(time, my_model.addPathogensToHosts, ["isolatedGV", { pos_al[-1]*5: 2, pos_al[-2]*5: 2}])

my_model.run(0,200) # Run the simulation for 200 time units.

data = my_model.saveToDataFrame('IsolatedIncreasingGV.csv')
    # Save the model results to a table.
graph13 = my_model.compartmentPlot('IsolatedIncreasingGV_HostsLT.png', data, populations=["isolatedGV"])
    # Plot the number of susceptible and infected hosts in the model over time.
graph21 = my_model.compositionPlot('ConnectedGV_GenotypesLT.png', data, populations=["population_mediumGV", "population_highGV"])
graph23 = my_model.compositionPlot('IsolatedIncreasingGV_GenotypesLT.png',data, populations=["isolatedGV"],
                                   track_specific_sequences=['MMMMM'])
    #Plot the track genotypes across time
graph3= my_model.populationsPlot('IsolatedIncreasingGVr_PopulationLT',data)
    #Track the number of infected, naive, death and recovered per population over time