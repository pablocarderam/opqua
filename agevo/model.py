#TODO: add non vector-borne mode, move gillespie and save model calls to gillespie.py

import numpy as np
import pandas as pd
from model import *
from gillespie import *

class Model(object):
    """docstring for Model."""

    def __init__(self):
        super(Model, self).__init__()
        self.populations = []
        self.setups = {}

    def newSetup(name,
        genome_length, genome_alphabet, treatments, contact_rate, treatment_rate,
        num_hosts, num_vectors, fitnessHost, fitnessVector,
        inoculum_host, inoculum_vector, inoculation_rate_host, inoculation_rate_vector,
        recovery_rate_host, recovery_rate_vector):
        """ Creates a new Parameters object with model parameters, saves it in
            setups dict under given name """

        self.setups[name] = Parameters(
            genome_length, genome_alphabet, treatments, contact_rate, treatment_rate,
            num_hosts, num_vectors, fitnessHost, fitnessVector,
            inoculum_host, inoculum_vector, inoculation_rate_host, inoculation_rate_vector,
            recovery_rate_host, recovery_rate_vector)


    def createPopulation(setup_name):
        """ Creates a new Population object with model parameters and adds it to
            the model. """

        self.populations.append( Population(self,self.setups[setup_name]) )


    def linkPopulations(pop1_index,pop2_index,rate):
        """ Establishes the migration rate between populations 1 and 2. """

        self.populations[pop1_index].neighbors[ self.populations[pop2_index] ] = rate


    def createInterconnectedPopulations(num_populations, setup_name, migration_rate):
        """ Creates new Population objects with model parameters and adds them to
            the model; links all of them to each other with same migration rate. """


        new_pops = [ Population(self,self.setups[setup_name]) for _ in range(num_populations) ]
        for i in range( len(self.populations), len(self.populations)+num_populations ):
            for j in range( len(self.populations), len(self.populations)+num_populations ):
                self.linkPopulations(i,j,rate)


        self.populations = self.populations + new_pops



    def run(time,save_to_dir=""):
        """ Runs the model until the given time, returns pandas dataframe with
            all data and saves it to a file if filepath given. """

        data = pd.DataFrame()

        def saveModelState():
            """ Saves model's current state to the pandas dataframe given,
                returns the modified dataframe. """

            return data

        def gillespieStep(arg):
            """ Calculates and advances model to the next Gillespie time step. """

            return pops



        if len(save_to_dir) > 0:
            data.write_csv(save_to_dir)

        return data
