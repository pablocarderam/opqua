
import numpy as np
import pandas as pd
from agevo.classes import *
from agevo.gillespie import *

class Model(object):
    """docstring for Model."""

    def __init__(self):
        super(Model, self).__init__()
        self.populations = {}
        self.setups = {}
        self.interventions = []

    def newSetup(name,
        num_loci, alleles_per_locus,
        num_hosts, num_vectors, fitnessHost, fitnessVector,
        contact_rate_host_vector, contact_rate_host_host,
        inoculum_host, inoculum_vector, inoculation_rate_host, inoculation_rate_vector,
        recovery_rate_host, recovery_rate_vector,
        recombine_in_host=0, recombine_in_vector=1,
        mutate_in_host=1e-6, mutate_in_vector=0,
        vector_borne=True, host_host_transmission=False):
        """ Creates a new Parameters object with model parameters, saves it in
            setups dict under given name """

        self.setups[name] = Parameters(
            num_loci, alleles_per_locus,
            num_hosts, num_vectors, fitnessHost, fitnessVector,
            contact_rate_host_vector, contact_rate_host_host,
            inoculum_host, inoculum_vector, inoculation_rate_host, inoculation_rate_vector,
            recovery_rate_host, recovery_rate_vector,
            recombine_in_host, recombine_in_vector,
            mutate_in_host, mutate_in_vector,
            vector_borne, host_host_transmission)


    def newIntervention(time, function, args):
        """ Creates a new Parameters object with model parameters, saves it in
            setups dict under given name """

        self.interventions.append( Intervention(time, function, args) )


    def run(t0,tf,save_to_dir=""):
        """ Runs the model between the given times, returns pandas dataframe with
            all data and saves it to a file if filepath given. """

        sim = Gillespie(self)
        data = sim.run(t0,tf)

        if len(save_to_dir) > 0:
            data.write_csv(save_to_dir)

        return data

    def createPopulation(setup_name,id):
        """ Creates a new Population object with model parameters and adds it to
            the model. """

        if id in self.populations.keys():
            id = id+'_2'

        self.populations[id] = Population(self,self.setups[setup_name],id)


    def linkPopulations(pop1_id,pop2_id,rate):
        """ Establishes the migration rate between populations 1 and 2. """

        self.populations[pop1_id].neighbors[ self.populations[pop2_id] ] = rate


    def createInterconnectedPopulations(num_populations, id_prefix, setup_name, migration_rate):
        """ Creates new Population objects with model parameters and adds them to
            the model; links all of them to each other with same migration rate. """


        new_pops = [ Population( self, self.setups[setup_name], id_prefix + str(i) ) for i in range(num_populations) ]
        for p1 in new_pops:
            if p1.id in self.populations.keys():
                p1.id = p1.id+'_2'

            self.populations[p1.id] = p1
            for p2 in new_pops:
                self.linkPopulations(p1.id,p2.id,rate)




    def addHosts(pop_id, num_hosts):
        """ Add a number of healthy hosts to population. """

        self.populations[pop_id].addHosts(num_hosts)


    def addVectors(pop_id, num_vectors):
        """ Add a number of healthy vectors to population """

        self.populations[pop_id].addVectors(num_vectors)


    def removeHosts(pop_id, num_hosts):
        """ Remove a number of random hosts from population """

        self.populations[pop_id].removeHosts(num_hosts)


    def removeVectors(pop_id, num_vectors):
        """ Remove a number of random vectors from population """

        self.populations[pop_id].removeVectors(num_vectors)


    def addPathogens( pop_id, strains, hosts=True ):
        """ Seeds pathogens according to strains dict (keys=genomes,
            values=num of infections); seeds on hosts unless hosts=False """

        self.populations[pop_id].addPathogens(strains,hosts)

    def treatHosts(pop_id, frac_hosts, treatment_seqs):
        """ Treat random hosts """

        self.populations[pop_id].treatHosts(frac_hosts,treatment_seqs)

    def treatVectors(pop_id, frac_vectors, treatment_seqs):
        """ Treat random vectors """

        self.populations[pop_id].treatVectors(frac_vectors,treatment_seqs)

    def protectHosts(pop_id, frac_hosts, protection_sequence):
        """ Treat random hosts """

        self.populations[pop_id].protectHosts(frac_hosts,protection_sequence)

    def protectVectors(pop_id, frac_vectors, protection_sequence):
        """ Treat random vectors """

        self.populations[pop_id].protectVectors(frac_vectors,protection_sequence)
