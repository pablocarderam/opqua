#TODO: select hosts/vectors and intervene on selections

import numpy as np
import pandas as pd
import difflib as dl
from opqua.internal.classes import *
from opqua.internal.gillespie import *

class Model(object):
    """docstring for Model."""

    def __init__(self):
        super(Model, self).__init__()
        self.populations = {}
        self.setups = {}
        self.interventions = []

    def newSetup(self, name, default=None,
        num_loci=None, possible_alleles=None,
        num_hosts=None, num_vectors=None, fitnessHost=None, fitnessVector=None,
        contact_rate_host_vector=None, contact_rate_host_host=None,
        inoculum_host=None, inoculum_vector=None, inoculation_rate_host=None, inoculation_rate_vector=None,
        recovery_rate_host=None, recovery_rate_vector=None,
        recombine_in_host=None, recombine_in_vector=None,
        mutate_in_host=None, mutate_in_vector=None):
        """ Creates a new Setup object with model parameters, saves it in
            setups dict under given name """

        if default == "vector-borne":
            num_loci = num_loci or 10
            possible_alleles = possible_alleles or 'ATCG'
            num_hosts = num_hosts or 100
            num_vectors = num_vectors or 100
            fitnessHost = fitnessHost or (lambda g: 1)
            fitnessVector = fitnessVector or (lambda g: 1)
            contact_rate_host_vector = contact_rate_host_vector or 1e-2
            contact_rate_host_host = contact_rate_host_host or 0
            inoculum_host = inoculum_host or 1e2
            inoculum_vector = inoculum_vector or 1e2
            inoculation_rate_host = inoculation_rate_host or 1e-1
            inoculation_rate_vector = inoculation_rate_vector or 1e-1
            recovery_rate_host = recovery_rate_host or 1e1
            recovery_rate_vector = recovery_rate_vector or 1e1
            recombine_in_host = recombine_in_host or 0
            recombine_in_vector = recombine_in_vector or 1e-2
            mutate_in_host = mutate_in_host or 1e-6
            mutate_in_vector = mutate_in_vector or 0

        elif default == "host-host":
            num_loci = num_loci or 10
            possible_alleles = possible_alleles or 'ATCG'
            num_hosts = num_hosts or 100
            num_vectors = num_vectors or 100
            fitnessHost = fitnessHost or (lambda g: 1)
            fitnessVector = fitnessVector or (lambda g: 1)
            contact_rate_host_vector = contact_rate_host_vector or 0
            contact_rate_host_host = contact_rate_host_host or 1e-2
            inoculum_host = inoculum_host or 1e2
            inoculum_vector = inoculum_vector or 1e2
            inoculation_rate_host = inoculation_rate_host or 1e-1
            inoculation_rate_vector = inoculation_rate_vector or 1e-1
            recovery_rate_host = recovery_rate_host or 5e-1
            recovery_rate_vector = recovery_rate_vector or 1e1
            recombine_in_host = recombine_in_host or 1e-3
            recombine_in_vector = recombine_in_vector or 0
            mutate_in_host = mutate_in_host or 1e-6
            mutate_in_vector = mutate_in_vector or 0

        self.setups[name] = Setup(
            num_loci, possible_alleles,
            num_hosts, num_vectors, fitnessHost, fitnessVector,
            contact_rate_host_vector, contact_rate_host_host,
            inoculum_host, inoculum_vector, inoculation_rate_host, inoculation_rate_vector,
            recovery_rate_host, recovery_rate_vector,
            recombine_in_host, recombine_in_vector,
            mutate_in_host, mutate_in_vector)


    def newIntervention(self, time, function, args):
        """ Creates a new Setup object with model parameters, saves it in
            setups dict under given name """

        self.interventions.append( Intervention(time, function, args) )


    def run(self,t0,tf,save_to_file):
        """ Runs the model between the given times, returns pandas dataframe with
            all data and saves it to a file if filepath given. """

        sim = Gillespie(self)
        data = sim.run(t0,tf,save_to_file)

        return data

    def newPopulation(self, id, setup_name, num_hosts=100, num_vectors=100):
        """ Creates a new Population object with model parameters and adds it to
            the model. """

        if id in self.populations.keys():
            id = id+'_2'

        self.populations[id] = Population( self, id, self.setups[setup_name], num_hosts, num_vectors )


    def linkPopulations(self, pop1_id, pop2_id, rate):
        """ Establishes the migration rate between populations 1 and 2. """

        self.populations[pop1_id].neighbors[ self.populations[pop2_id] ] = rate


    def createInterconnectedPopulations(self, num_populations, migration_rate, id_prefix, setup_name, num_hosts=100, num_vectors=100):
        """ Creates new Population objects with model parameters and adds them to
            the model; links all of them to each other with same migration rate. """


        new_pops = [ Population( self, id_prefix + str(i), self.setups[setup_name], num_hosts, num_vectors ) for i in range(num_populations) ]
        for p1 in new_pops:
            if p1.id in self.populations.keys():
                p1.id = p1.id+'_2'

            self.populations[p1.id] = p1
            for p2 in new_pops:
                self.linkPopulations(p1.id,p2.id,rate)




    def addHosts(self, pop_id, num_hosts):
        """ Add a number of healthy hosts to population. """

        return self.populations[pop_id].addHosts(num_hosts)


    def addVectors(self, pop_id, num_vectors):
        """ Add a number of healthy vectors to population """

        return self.populations[pop_id].addVectors(num_vectors)


    def removeHosts(self, pop_id, num_hosts_or_list):
        """ Remove a number of random hosts from population """

        self.populations[pop_id].removeHosts(num_hosts_or_list)


    def removeVectors(self, pop_id, num_vectors_or_list):
        """ Remove a number of random vectors from population """

        self.populations[pop_id].removeVectors(num_vectors_or_list)


    def addPathogensToHosts(self, pop_id, genomes_numbers, hosts=[]):
        """ Seeds pathogens according to strains dict (keys=genomes,
            values=num of infections); seeds on hosts unless hosts=False """

        self.populations[pop_id].addPathogensToHosts(genomes_numbers,hosts)

    def addPathogensToVectors(self, pop_id, genomes_numbers, vectors=[]):
        """ Seeds pathogens according to strains dict (keys=genomes,
            values=num of infections); seeds on hosts unless hosts=False """

        self.populations[pop_id].addPathogensToVectors(genomes_numbers,vectors)

    def treatHosts(self, pop_id, frac_hosts, treatment_seqs, hosts=[]):
        """ Treat random hosts """

        self.populations[pop_id].treatHosts(frac_hosts,treatment_seqs,hosts)

    def treatVectors(self, pop_id, frac_vectors, treatment_seqs, vectors=[]):
        """ Treat random vectors """

        self.populations[pop_id].treatVectors(frac_vectors,treatment_seqs,vectors)

    def protectHosts(self, pop_id, frac_hosts, protection_sequence, hosts=[]):
        """ Treat random hosts """

        self.populations[pop_id].protectHosts(frac_hosts,protection_sequence,hosts)

    def protectVectors(self, pop_id, frac_vectors, protection_sequence, vectors=[]):
        """ Treat random vectors """

        self.populations[pop_id].protectVectors(frac_vectors,protection_sequence,vectors)

    def setSetup(self, pop_id, setup_id):
        """ Treat random vectors """

        self.populations[pop_id].setSetup( self.setups[setup_id] )

    # Fitness functions
    @staticmethod
    def stabilizingSelection(genome,optimal_genome,min_fitness):
        """ A purifying selection fitness function based on exponential decay of
            fitness based on distance from an optimal genome sequence. """

        similarity = dl.SequenceMatcher(None, genome, optimal_genome).ratio()
        fitness = np.exp( np.log( min_fitness ) * ( 1-similarity ) )

        return fitness

    @staticmethod
    def disruptiveSelection(genome,worst_genome,min_fitness):
        """ A purifying selection fitness function based on exponential decay of
            fitness based on distance from an optimal genome sequence. """

        similarity = dl.SequenceMatcher(None, genome, optimal_genome).ratio()
        fitness = np.exp( np.log( min_fitness ) * ( similarity ) )

        return fitness
