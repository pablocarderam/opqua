
"""Contains class Model; main class user interacts with."""

import numpy as np
import pandas as pd
import textdistance as td
import itertools as it
import copy as cp
import seaborn as sns
import joblib as jl

from opqua.internal.host import Host
from opqua.internal.vector import Vector
from opqua.internal.population import Population
from opqua.internal.setup import Setup
from opqua.internal.intervention import Intervention
from opqua.internal.gillespie import Gillespie
from opqua.internal.data import saveToDf, getPathogens, getProtections, \
    getPathogenDistanceHistoryDf
from opqua.internal.plot import populationsPlot, compartmentPlot, \
    compositionPlot, clustermap

class Model(object):
    """Class defines a Model.

    This is the main class that the user interacts with.

    The Model class contains populations, setups, and interventions to be used
    in simulation. Also contains groups of hosts/vectors for manipulations and
    stores model history as snapshots for each time point.

    **CONSTANTS:**

    - `CB_PALETTE`: a colorblind-friendly 8-color color scheme.
    - `DEF_CMAP`: a colormap object for Seaborn plots.

    Attributes:
        populations: dictionary with keys=population IDs, values=Population
            objects.
        setups: dictionary with keys=setup IDs, values=Setup objects.
        interventions: contains model interventions in the order they will occur.
        groups: dictionary with keys=group IDs, values=lists of hosts/vectors.
        history: dictionary with keys=time values, values=Model objects that
            are snapshots of Model at that timepoint.
        t_var: variable that tracks time in simulations.
    """

    ### CONSTANTS ###
    ### Color scheme constants ###
    CB_PALETTE = ["#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]
     # www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
     # http://jfly.iam.u-tokyo.ac.jp/color/

    DEF_CMAP = sns.cubehelix_palette(
        start=.5, rot=-.75, as_cmap=True, reverse=True
        )


    ### CLASS CONSTRUCTOR ###

    def __init__(self):
        """Create a new Model object."""
        super(Model, self).__init__()
        self.populations = {}
            # dictionary with keys=population IDs, values=Population objects
        self.setups = {}
            # dictionary with keys=setup IDs, values=Setup objects
        self.interventions = []
            # contains model interventions in the order they will occur
        self.groups = {}
            # dictionary with keys=group IDs, values=lists of hosts/vectors
        self.history = {}
            # dictionary with keys=time values, values=Model objects that are
            # snapshots of Model at that timepoint
        self.global_trackers = {
                # dictionary keeping track of some global indicators over all
                # the course of the simulation
            'num_events' : { id:0 for id in Gillespie.EVENT_IDS.values() },
                # tracks the number of each kind of event in the simulation
            'last_event_time' : 0,
                # time point at which the last event in the simulation happened
            'genomes_seen' : [],
                # list of all unique genomes that have appeared in the
                # simulation
            'custom_conditions' : {}
                # dictionary with keys=ID of custom condition, values=lists of
                # times; every time True is returned by a function in
                # custom_condition_trackers, the simulation time will be stored
                # under the corresponding ID inside
                # global_trackers['custom_condition']
            }
        self.custom_condition_trackers = {}
            # dictionary with keys=ID of custom condition, values=functions that
            # take a Model object as argument and return True or False; every
            # time True is returned by a function in custom_condition_trackers,
            # the simulation time will be stored under the corresponding ID
            # inside global_trackers['custom_condition']

        self.t_var = 0 # used as time variable during simulation

    ### MODEL METHODS ###

    ### Model initialization and simulation: ###

    def setRandomSeed(self, seed):
        """Set random seed for numpy random number generator.

        Arguments:
            seed (int): int for the random seed to be passed to numpy.
        """

        np.random.seed(seed)

    def newSetup(
            self, name, preset=None,
            num_loci=None, possible_alleles=None,
            fitnessHost=None, contactHost=None, receiveContactHost=None,
            mortalityHost=None, natalityHost=None,
            recoveryHost=None, migrationHost=None,
            populationContactHost=None, receivePopulationContactHost=None,
            mutationHost=None,
            recombinationHost=None, fitnessVector=None,
            contactVector=None, receiveContactVector=None, mortalityVector=None,
            natalityVector=None, recoveryVector=None,
            migrationVector=None, populationContactVector=None,
            receivePopulationContactVector=None,
            mutationVector=None, recombinationVector=None,
            contact_rate_host_vector=None,
            transmission_efficiency_host_vector=None,
            transmission_efficiency_vector_host=None,
            contact_rate_host_host=None,
            transmission_efficiency_host_host=None,
            mean_inoculum_host=None, mean_inoculum_vector=None,
            recovery_rate_host=None, recovery_rate_vector=None,
            mortality_rate_host=None, mortality_rate_vector=None,
            recombine_in_host=None, recombine_in_vector=None,
            num_crossover_host=None, num_crossover_vector=None,
            mutate_in_host=None, mutate_in_vector=None,
            death_rate_host=None, death_rate_vector=None,
            birth_rate_host=None, birth_rate_vector=None,
            vertical_transmission_host=None, vertical_transmission_vector=None,
            inherit_protection_host=None, inherit_protection_vector=None,
            protection_upon_recovery_host=None,
            protection_upon_recovery_vector=None):
        """Create a new `Setup`, save it in setups dict under given name.

        Two preset setups exist: "vector-borne" and "host-host". You may select
        one of the preset setups with the preset keyword argument and then
        modify individual parameters with additional keyword arguments, without
        having to specify all of them.

        **"host-host":**

        - `num_loci` = 10
        - `possible_alleles` = 'ATCG'
        - `fitnessHost` = (lambda g: 1)
        - `contactHost` = (lambda g: 1)
        - `receiveContactHost` = (lambda g: 1)
        - `mortalityHost` = (lambda g: 1)
        - `natalityHost` = (lambda g: 1)
        - `recoveryHost` = (lambda g: 1)
        - `migrationHost` = (lambda g: 1)
        - `populationContactHost` = (lambda g: 1)
        - `receivePopulationContactHost` = (lambda g: 1)
        - `mutationHost` = (lambda g: 1)
        - `recombinationHost` = (lambda g: 1)
        - `fitnessVector` = (lambda g: 1)
        - `contactVector` = (lambda g: 1)
        - `receiveContactVector` = (lambda g: 1)
        - `mortalityVector` = (lambda g: 1)
        - `natalityVector` = (lambda g: 1)
        - `recoveryVector` = (lambda g: 1)
        - `migrationVector` = (lambda g: 1)
        - `populationContactVector` = (lambda g: 1)
        - `receivePopulationContactVector` = (lambda g: 1)
        - `mutationVector` = (lambda g: 1)
        - `recombinationVector` = (lambda g: 1)
        - `contact_rate_host_vector` = 0
        - `transmission_efficiency_host_vector` = 0
        - `transmission_efficiency_vector_host` = 0
        - `contact_rate_host_host` = 2e-1
        - `transmission_efficiency_host_host` = 1
        - `mean_inoculum_host` = 1e1
        - `mean_inoculum_vector` = 0
        - `recovery_rate_host` = 1e-1
        - `recovery_rate_vector` = 0
        - `mortality_rate_host` = 0
        - `mortality_rate_vector` = 0
        - `recombine_in_host` = 1e-4
        - `recombine_in_vector` = 0
        - `num_crossover_host` = 1
        - `num_crossover_vector` = 0
        - `mutate_in_host` = 1e-6
        - `mutate_in_vector` = 0
        - `death_rate_host` = 0
        - `death_rate_vector` = 0
        - `birth_rate_host` = 0
        - `birth_rate_vector` = 0
        - `vertical_transmission_host` = 0
        - `vertical_transmission_vector` = 0
        - `inherit_protection_host` = 0
        - `inherit_protection_vector` = 0
        - `protection_upon_recovery_host` = None
        - `protection_upon_recovery_vector` = None

        **"vector-borne":**

        - `num_loci` = 10
        - `possible_alleles` = 'ATCG'
        - `fitnessHost` = (lambda g: 1)
        - `contactHost` = (lambda g: 1)
        - `receiveContactHost` = (lambda g: 1)
        - `mortalityHost` = (lambda g: 1)
        - `natalityHost` = (lambda g: 1)
        - `recoveryHost` = (lambda g: 1)
        - `migrationHost` = (lambda g: 1)
        - `populationContactHost` = (lambda g: 1)
        - `receivePopulationContactHost` = (lambda g: 1)
        - `mutationHost` = (lambda g: 1)
        - `recombinationHost` = (lambda g: 1)
        - `fitnessVector` = (lambda g: 1)
        - `contactVector` = (lambda g: 1)
        - `receiveContactVector` = (lambda g: 1)
        - `mortalityVector` = (lambda g: 1)
        - `natalityVector` = (lambda g: 1)
        - `recoveryVector` = (lambda g: 1)
        - `migrationVector` = (lambda g: 1)
        - `populationContactVector` = (lambda g: 1)
        - `receivePopulationContactVector` = (lambda g: 1)
        - `mutationVector` = (lambda g: 1)
        - `recombinationVector` = (lambda g: 1)
        - `contact_rate_host_vector` = 2e-1
        - `transmission_efficiency_host_vector` = 1
        - `transmission_efficiency_vector_host` = 1
        - `contact_rate_host_host` = 0
        - `transmission_efficiency_host_host` = 0
        - `mean_inoculum_host` = 1e2
        - `mean_inoculum_vector` = 1e0
        - `recovery_rate_host` = 1e-1
        - `recovery_rate_vector` = 1e-1
        - `mortality_rate_host` = 0
        - `mortality_rate_vector` = 0
        - `recombine_in_host` = 0
        - `recombine_in_vector` = 1e-4
        - `num_crossover_host` = 0
        - `num_crossover_vector` = 1
        - `mutate_in_host` = 1e-6
        - `mutate_in_vector` = 0
        - `death_rate_host` = 0
        - `death_rate_vector` = 0
        - `birth_rate_host` = 0
        - `birth_rate_vector` = 0
        - `vertical_transmission_host` = 0
        - `vertical_transmission_vector` = 0
        - `inherit_protection_host` = 0
        - `inherit_protection_vector` = 0
        - `protection_upon_recovery_host` = None
        - `protection_upon_recovery_vector` = None

        Arguments:
            name (String): name of setup to be used as a key in model setups dictionary.

        Keyword arguments:
            preset (None or String): preset setup to be used: "vector-borne" or "host-host", if
                None, must define all other keyword arguments. Defaults to None.
            num_loci (int>0): length of each pathogen genome string.
            possible_alleles (String or list of Strings with num_loci elements): set of possible 
                characters in all genome string, or at each position in genome string. 
            fitnessHost (callable, takes a String argument and returns a number >= 0): function 
                that evaluates relative fitness in head-to-head competition for different genomes 
                within the same host.
            contactHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying probability of a given host being chosen to be the 
                infector in a contact event, based on genome sequence of pathogen.
            receiveContactHost (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying probability of a given host being chosen to be 
                the infected in a contact event, based on genome sequence of pathogen.
            mortalityHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying death rate for a given host, based on genome sequence 
                of pathogen.
            natalityHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying birth rate for a given host, based on genome sequence 
                of pathogen.
            recoveryHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying recovery rate for a given host based on genome sequence 
                of pathogen.
            migrationHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying migration rate for a given host based on genome sequence 
                of pathogen.
            populationContactHost (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying population contact rate for a given host based on 
                genome sequence of pathogen.
            mutationHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying mutation rate for a given host based on genome sequence 
                of pathogen.
            recombinationHost (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying recombination rate for a given host based on genome 
                sequence of pathogen.
            fitnessVector (callable, takes a String argument and returns a number >=0): function that 
                evaluates relative fitness in head-to-head competition for different genomes within 
                the same vector.
            contactVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying probability of a given vector being chosen to be the 
                infector in a contact event, based on genome sequence of pathogen.
            receiveContactVector (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying probability of a given vector being chosen to be the 
                infected in a contact event, based on genome sequence of pathogen.
            mortalityVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying death rate for a given vector, based on genome sequence 
                of pathogen.
            natalityVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying birth rate for a given vector, based on genome sequence 
                of pathogen.
            recoveryVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying recovery rate for a given vector based on genome sequence 
                of pathogen.
            migrationVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying migration rate for a given vector based on genome sequence 
                of pathogen.
            populationContactVector (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying population contact rate for a given vector based on 
                genome sequence of pathogen.
            mutationVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying mutation rate for a given vector based on genome sequence 
                of pathogen.
            recombinationVector (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying recombination rate for a given vector based on genome 
                sequence of pathogen.
            contact_rate_host_vector (number >= 0): rate of host-vector contact events, not necessarily 
                transmission, assumes constant population density; evts/time.
            transmission_efficiency_host_vector (float): fraction of host-vector contacts
                that result in successful transmission.
            transmission_efficiency_vector_host (float): fraction of vector-host contacts
                that result in successful transmission.
            contact_rate_host_host (number >= 0): rate of host-host contact events, not
                necessarily transmission, assumes constant population density; evts/time.
            transmission_efficiency_host_host (float): fraction of host-host contacts
                    that result in successful transmission.
            mean_inoculum_host (int >= 0): mean number of pathogens that are transmitted from
                a vector or host into a new host during a contact event.
            mean_inoculum_vector (int >= 0) mean number of pathogens that are transmitted
                from a host to a vector during a contact event.
            recovery_rate_host (number >= 0): rate at which hosts clear all pathogens;
                1/time.
            recovery_rate_vector (number >= 0): rate at which vectors clear all pathogens
                1/time.
            recovery_rate_vector (number >= 0): rate at which vectors clear all pathogens
                1/time.
            mortality_rate_host (number 0-1): rate at which infected hosts die from disease.
            mortality_rate_vector (number 0-1): rate at which infected vectors die from
                disease.
            recombine_in_host (number >= 0): rate at which recombination occurs in host;
                evts/time.
            recombine_in_vector (number >= 0): rate at which recombination occurs in vector;
                evts/time.
            num_crossover_host (number >= 0): mean of a Poisson distribution modeling the number
                of crossover events of host recombination events.
            num_crossover_vector (number >= 0): mean of a Poisson distribution modeling the
                number of crossover events of vector recombination events.
            mutate_in_host (number >= 0): rate at which mutation occurs in host; evts/time.
            mutate_in_vector (number >= 0): rate at which mutation occurs in vector; evts/time.
            death_rate_host (number >= 0): natural host death rate; 1/time.
            death_rate_vector (number >= 0): natural vector death rate; 1/time.
            birth_rate_host (number >= 0): infected host birth rate; 1/time.
            birth_rate_vector (number >= 0): infected vector birth rate; 1/time.
            vertical_transmission_host (number 0-1): probability that a host is infected by its
                parent at birth.
            vertical_transmission_vector (number 0-1): probability that a vector is infected by
                its parent at birth.
            inherit_protection_host (number 0-1): probability that a host inherits all
                protection sequences from its parent.
            inherit_protection_vector (number 0-1): probability that a vector inherits all
                protection sequences from its parent.
            protection_upon_recovery_host (None or array-like of length 2 with int 0-num_loci): defines 
                indexes in genome string that define substring to be added to host protection sequences 
                after recovery.
            protection_upon_recovery_vector (None or array-like of length 2 with int 0-num_loci): defines 
                indexes in genome string that define substring to be added to vector protection sequences 
                after recovery.
        """

        if preset == "vector-borne":
            num_loci = 10 if num_loci is None else num_loci
            possible_alleles = \
                'ATCG' if possible_alleles is None else possible_alleles
            fitnessHost = (lambda g: 1) if fitnessHost is None else fitnessHost
            contactHost = (lambda g: 1) if contactHost is None else contactHost
            receiveContactHost = \
                (lambda g: 1) if receiveContactHost is None else receiveContactHost
            mortalityHost = \
                (lambda g: 1) if mortalityHost is None else mortalityHost
            natalityHost = \
                (lambda g: 1) if natalityHost is None else natalityHost
            recoveryHost = \
                (lambda g: 1) if recoveryHost is None else recoveryHost
            migrationHost = \
                (lambda g: 1) if migrationHost is None else migrationHost
            populationContactHost = \
                (lambda g: 1) if populationContactHost is None else populationContactHost
            receivePopulationContactHost = \
                (lambda g: 1) if receivePopulationContactHost is None else receivePopulationContactHost
            mutationHost = \
                (lambda g: 1) if mutationHost is None else mutationHost
            recombinationHost = \
                (lambda g: 1) if recombinationHost is None else recombinationHost
            fitnessVector = \
                (lambda g: 1) if fitnessVector is None else fitnessVector
            contactVector = \
                (lambda g: 1) if contactVector is None else contactVector
            receiveContactVector = \
                (lambda g: 1) if receiveContactVector is None else receiveContactVector
            mortalityVector = \
                (lambda g: 1) if mortalityVector is None else mortalityVector
            natalityVector = \
                (lambda g: 1) if natalityVector is None else natalityVector
            recoveryVector = \
                (lambda g: 1) if recoveryVector is None else recoveryVector
            migrationVector = \
                (lambda g: 1) if migrationVector is None else migrationVector
            populationContactVector = \
                (lambda g: 1) if populationContactVector is None else populationContactVector
            receivePopulationContactVector = \
                (lambda g: 1) if receivePopulationContactVector is None else receivePopulationContactVector
            mutationVector = \
                (lambda g: 1) if mutationVector is None else mutationVector
            recombinationVector = \
                (lambda g: 1) if recombinationVector is None else recombinationVector
            contact_rate_host_vector = \
                2e-1 if contact_rate_host_vector is None \
                else contact_rate_host_vector
            transmission_efficiency_host_vector = \
                1 if transmission_efficiency_host_vector is None \
                else transmission_efficiency_host_vector
            transmission_efficiency_vector_host = \
                1 if transmission_efficiency_vector_host is None \
                else transmission_efficiency_vector_host
            contact_rate_host_host = \
                0 if contact_rate_host_host is None else contact_rate_host_host
            transmission_efficiency_host_host = \
                0 if transmission_efficiency_host_host is None \
                else transmission_efficiency_host_host
            mean_inoculum_host = \
                1e2 if mean_inoculum_host is None else mean_inoculum_host
            mean_inoculum_vector = \
                1 if mean_inoculum_vector is None else mean_inoculum_vector
            recovery_rate_host = \
                1e-1 if recovery_rate_host is None else recovery_rate_host
            recovery_rate_vector = \
                1e-1 if recovery_rate_vector is None else recovery_rate_vector
            mortality_rate_host = \
                0 if mortality_rate_host is None else mortality_rate_host
            mortality_rate_vector = \
                0 if mortality_rate_vector is None else mortality_rate_vector
            recombine_in_host = \
                0 if recombine_in_host is None else recombine_in_host
            recombine_in_vector = \
                1e-4 if recombine_in_vector is None else recombine_in_vector
            num_crossover_host = 0 \
                if num_crossover_host is None else num_crossover_host
            num_crossover_vector = \
                1 if num_crossover_vector is None else num_crossover_vector
            mutate_in_host = 1e-6 if mutate_in_host is None else mutate_in_host
            mutate_in_vector = \
                0 if mutate_in_vector is None else mutate_in_vector
            death_rate_host = 0 if death_rate_host is None else death_rate_host
            death_rate_vector = \
                0 if death_rate_vector is None else death_rate_vector
            birth_rate_host = 0 if birth_rate_host is None else birth_rate_host
            birth_rate_vector = \
                0 if birth_rate_vector is None else birth_rate_vector
            vertical_transmission_host = \
                0 if vertical_transmission_host is None \
                else vertical_transmission_host
            vertical_transmission_vector = \
                0 if vertical_transmission_vector is None \
                else vertical_transmission_vector
            inherit_protection_host = \
                0 if inherit_protection_host is None \
                else inherit_protection_host
            inherit_protection_vector = \
                0 if inherit_protection_vector is None \
                else inherit_protection_vector
            protection_upon_recovery_host = protection_upon_recovery_host
            protection_upon_recovery_vector = protection_upon_recovery_vector

        elif preset == "host-host":
            num_loci = 10 if num_loci is None else num_loci
            possible_alleles = \
                'ATCG' if possible_alleles is None else possible_alleles
            fitnessHost = (lambda g: 1) if fitnessHost is None else fitnessHost
            contactHost = (lambda g: 1) if contactHost is None else contactHost
            receiveContactHost = \
                (lambda g: 1) if receiveContactHost is None else receiveContactHost
            mortalityHost = \
                (lambda g: 1) if mortalityHost is None else mortalityHost
            natalityHost = \
                (lambda g: 1) if natalityHost is None else natalityHost
            recoveryHost = \
                (lambda g: 1) if recoveryHost is None else recoveryHost
            migrationHost = \
                (lambda g: 1) if migrationHost is None else migrationHost
            populationContactHost = \
                (lambda g: 1) if populationContactHost is None else populationContactHost
            receivePopulationContactHost = \
                (lambda g: 1) if receivePopulationContactHost is None else receivePopulationContactHost
            mutationHost = \
                (lambda g: 1) if mutationHost is None else mutationHost
            recombinationHost = \
                (lambda g: 1) if recombinationHost is None else recombinationHost
            fitnessVector = \
                (lambda g: 1) if fitnessVector is None else fitnessVector
            contactVector = \
                (lambda g: 1) if contactVector is None else contactVector
            receiveContactVector = \
                (lambda g: 1) if receiveContactVector is None else receiveContactVector
            mortalityVector = \
                (lambda g: 1) if mortalityVector is None else mortalityVector
            natalityVector = \
                (lambda g: 1) if natalityVector is None else natalityVector
            recoveryVector = \
                (lambda g: 1) if recoveryVector is None else recoveryVector
            mortality_rate_host = \
                0 if mortality_rate_host is None else mortality_rate_host
            mortality_rate_vector = \
                0 if mortality_rate_vector is None else mortality_rate_vector
            migrationVector = \
                (lambda g: 1) if migrationVector is None else migrationVector
            populationContactVector = \
                (lambda g: 1) if populationContactVector is None else populationContactVector
            receivePopulationContactVector = \
                (lambda g: 1) if receivePopulationContactVector is None else receivePopulationContactVector
            mutationVector = \
                (lambda g: 1) if mutationVector is None else mutationVector
            recombinationVector = \
                (lambda g: 1) if recombinationVector is None else recombinationVector
            contact_rate_host_vector = \
                0 if contact_rate_host_vector is None \
                else contact_rate_host_vector
            transmission_efficiency_host_vector = \
                0 if transmission_efficiency_host_vector is None \
                else transmission_efficiency_host_vector
            transmission_efficiency_vector_host = \
                0 if transmission_efficiency_vector_host is None \
                else transmission_efficiency_vector_host
            contact_rate_host_host = \
                2e-1 if contact_rate_host_host is None \
                else contact_rate_host_host
            transmission_efficiency_host_host = \
                1 if transmission_efficiency_host_host is None \
                else transmission_efficiency_host_host
            mean_inoculum_host = \
                1e1 if mean_inoculum_host is None else mean_inoculum_host
            mean_inoculum_vector = \
                0 if mean_inoculum_vector is None else mean_inoculum_vector
            recovery_rate_host = \
                1e-1 if recovery_rate_host is None else recovery_rate_host
            recovery_rate_vector = \
                0 if recovery_rate_vector is None else recovery_rate_vector
            recombine_in_host = \
                1e-4 if recombine_in_host is None else recombine_in_host
            recombine_in_vector = \
                0 if recombine_in_vector is None else recombine_in_vector
            num_crossover_host = 1 \
                if num_crossover_host is None else num_crossover_host
            num_crossover_vector = \
                0 if num_crossover_vector is None else num_crossover_vector
            mutate_in_host = \
                1e-6 if mutate_in_host is None else mutate_in_host
            mutate_in_vector = \
                0 if mutate_in_vector is None else mutate_in_vector
            death_rate_host = \
                0 if death_rate_host is None else death_rate_host
            death_rate_vector = \
                0 if death_rate_vector is None else death_rate_vector
            birth_rate_host = 0 if birth_rate_host is None else birth_rate_host
            birth_rate_vector = \
                0 if birth_rate_vector is None else birth_rate_vector
            vertical_transmission_host = \
                0 if vertical_transmission_host is None \
                else vertical_transmission_host
            vertical_transmission_vector = \
                0 if vertical_transmission_vector is None \
                else vertical_transmission_vector
            inherit_protection_host = \
                0 if inherit_protection_host is None \
                else inherit_protection_host
            inherit_protection_vector = \
                0 if inherit_protection_vector is None \
                else inherit_protection_vector
            protection_upon_recovery_host = protection_upon_recovery_host
            protection_upon_recovery_vector = protection_upon_recovery_vector

        self.setups[name] = Setup(
            name,
            num_loci, possible_alleles,
            fitnessHost, contactHost, receiveContactHost, mortalityHost,
            natalityHost, recoveryHost, migrationHost,
            populationContactHost, receivePopulationContactHost,
            mutationHost, recombinationHost,
            fitnessVector, contactVector, receiveContactVector, mortalityVector,
            natalityVector,recoveryVector, migrationVector,
            populationContactVector, receivePopulationContactVector,
            mutationVector, recombinationVector,
            contact_rate_host_vector,
            transmission_efficiency_host_vector,
            transmission_efficiency_vector_host,
            contact_rate_host_host,
            transmission_efficiency_host_host,
            mean_inoculum_host, mean_inoculum_vector,
            recovery_rate_host, recovery_rate_vector,
            mortality_rate_host,mortality_rate_vector,
            recombine_in_host, recombine_in_vector,
            num_crossover_host, num_crossover_vector,
            mutate_in_host, mutate_in_vector, death_rate_host,death_rate_vector,
            birth_rate_host, birth_rate_vector,
            vertical_transmission_host, vertical_transmission_vector,
            inherit_protection_host, inherit_protection_vector,
            protection_upon_recovery_host, protection_upon_recovery_vector
            )

    def newIntervention(self, time, method_name, args):
        """Create a new intervention to be carried out at a specific time.

        Arguments:
            time (number >= 0): time at which intervention will take place.
            method_name (String): intervention to be carried out, must correspond to the
                name of a method of the `Model` object.
            args (array-like): contains arguments for function in positinal order.
        """

        self.interventions.append( Intervention(time, method_name, args, self) )

    def addCustomConditionTracker(self, condition_id, trackerFunction):
        """Add a function to track occurrences of custom events in simulation.

        Adds function `trackerFunction` to dictionary `custom_condition_trackers`
        under key `condition_id`. Function `trackerFunction` will be executed at
        every event in the simulation. Every time True is returned,
        the simulation time will be stored under the corresponding `condition_id`
        key inside `global_trackers['custom_condition']`.

        Arguments:
            condition_id (String): ID of this specific condition-
            trackerFunction (callable): function that take a `Model` object as argument 
                and returns True or False.
        """

        self.custom_condition_trackers['condition_id'] = trackerFunction
        self.global_trackers['custom_conditions']['condition_id'] = []

    def run(self,t0,tf,time_sampling=0,host_sampling=0,vector_sampling=0):
        """Simulate model for a specified time between two time points.

        Simulates a time series using the Gillespie algorithm.

        Saves a dictionary containing model state history, with `keys=times` and
        `values=Model` objects with model snapshot at that time point under this
        model's history attribute.

        Arguments:
            t0 (number >= 0): initial time point to start simulation at.
            tf (number >= 0): initial time point to end simulation at.

        Keyword arguments:
            time_sampling (int): how many events to skip before saving a snapshot of the
                system state (saves all by default), if <0, saves only final state. Defaults to 0.
            host_sampling (int >= 0): how many hosts to skip before saving one in a snapshot
                of the system state (saves all by default). Defaults to 0.
            vector_sampling (int >= 0): how many vectors to skip before saving one in a
                snapshot of the system state (saves all by default). Defaults to 0.
        """

        sim = Gillespie(self)
        self.history = sim.run(
            t0, tf, time_sampling, host_sampling, vector_sampling
            )

    def runReplicates(
            self,t0,tf,replicates,host_sampling=0,vector_sampling=0,n_cores=0,
            **kwargs):
        """Simulate replicates of a model, save only end results.

        Simulates replicates of a time series using the Gillespie algorithm.

        Saves a dictionary containing model end state state, with `keys=times` and
        `values=Model` objects with model snapshot. The time is the final timepoint.

        Arguments:
            t0 (number >= 0): initial time point to start simulation at.
            tf (number >= 0): initial time point to end simulation at.
            replicates (int >= 1): how many replicates to simulate.

        Keyword arguments:
            host_sampling (int >= 0): how many hosts to skip before saving one in a snapshot
                of the system state (saves all by default). Defaults to 0.
            vector_sampling (int >= 0): how many vectors to skip before saving one in a
                snapshot of the system state (saves all by default). Defaults to 0.
            n_cores (int >= 0): number of cores to parallelize file export across, if 0, all
                cores available are used. Defaults to 0.
            **kwargs: additional arguents for joblib multiprocessing.

        Returns:
            List of `Model` objects with the final snapshots.
        """

        if not n_cores:
            n_cores = jl.cpu_count()

        print('Starting parallel simulations...')

        def run(sim_num):
            model = self.deepCopy()
            sim = Gillespie(model)
            mod = sim.run(
                t0,tf,time_sampling=-1,
                host_sampling=host_sampling,vector_sampling=vector_sampling
                )[tf]
            mod.history = { tf:mod }
            return mod

        return jl.Parallel(n_jobs=n_cores, verbose=10, **kwargs) (
            jl.delayed( run ) (_) for _ in range(replicates)
            )

    def runParamSweep(
            self,t0,tf,setup_id,
            param_sweep_dic={},
            host_population_size_sweep={}, vector_population_size_sweep={},
            host_migration_sweep_dic={}, vector_migration_sweep_dic={},
            host_host_population_contact_sweep_dic={},
            host_vector_population_contact_sweep_dic={},
            vector_host_population_contact_sweep_dic={},
            replicates=1,host_sampling=0,vector_sampling=0,n_cores=0,
            **kwargs):
        """Simulate a parameter sweep with a model, save only end results.

        Simulates variations of a time series using the Gillespie algorithm.

        Saves a dictionary containing model end state state, with `keys=times` and
        `values=Model` objects with model snapshot. The time is the final
        timepoint.

        Arguments:
            t0 (number >= 0): initial time point to start simulation at.
            tf (number >= 0): initial time point to end simulation at.
            setup_id (String): ID of setup to be assigned.

        Keyword arguments:
            param_sweep_dic -- dictionary with keys=parameter names (attributes of
                Setup), values=list of values for parameter (list, class of elements
                depends on parameter)
            host_population_size_sweep -- dictionary with keys=population IDs
                (Strings), values=list of values with host population sizes
                (must be greater than original size set for each population, list of
                numbers)
            vector_population_size_sweep -- dictionary with keys=population IDs
                (Strings), values=list of values with vector population sizes
                (must be greater than original size set for each population, list of
                numbers)
            host_migration_sweep_dic -- dictionary with keys=population IDs of
                origin and destination, separated by a colon ';' (Strings),
                values=list of values (list of numbers)
            vector_migration_sweep_dic -- dictionary with keys=population IDs of
                origin and destination, separated by a colon ';' (Strings),
                values=list of values (list of numbers)
            host_host_population_contact_sweep_dic -- dictionary with
                keys=population IDs of origin and destination, separated by a colon
                ';' (Strings), values=list of values (list of numbers)
            host_vector_population_contact_sweep_dic -- dictionary with
                keys=population IDs of origin and destination, separated by a colon
                ';' (Strings), values=list of values (list of numbers)
            vector_host_population_contact_sweep_dic -- dictionary with
                keys=population IDs of origin and destination, separated by a colon
                ';' (Strings), values=list of values (list of numbers)
            replicates -- how many replicates to simulate (int >= 1)
            host_sampling -- how many hosts to skip before saving one in a snapshot
                of the system state (saves all by default) (int >= 0, default 0)
            vector_sampling -- how many vectors to skip before saving one in a
                snapshot of the system state (saves all by default)
                (int >= 0, default 0)
            n_cores -- number of cores to parallelize file export across, if 0, all
                cores available are used (default 0; int >= 0)
            **kwargs -- additional arguents for joblib multiprocessing

        Returns:
            DataFrame with parameter combinations, list of Model objects with the
                final snapshots.
        """

        if not n_cores:
            n_cores = jl.cpu_count()

        for p in host_population_size_sweep:
            param_sweep_dic['pop_size_host:'+p] = host_population_size_sweep[p]

        for p in vector_population_size_sweep:
            param_sweep_dic['pop_size_vector:'+p] = vector_population_size_sweep[p]

        for p in host_migration_sweep_dic:
            param_sweep_dic['migrate_host:'+p] = host_migration_sweep_dic[p]

        for p in vector_migration_sweep_dic:
            param_sweep_dic['migrate_vector:'+p] = vector_migration_sweep_dic[p]

        for p in host_host_population_contact_sweep_dic:
            param_sweep_dic['population_contact_host_host:'+p] \
                = host_host_population_contact_sweep_dic[p]

        for p in host_vector_population_contact_sweep_dic:
            param_sweep_dic['population_contact_host_vector:'+p] \
                = host_vector_population_contact_sweep_dic[p]

        for p in vector_host_population_contact_sweep_dic:
            param_sweep_dic['population_contact_vector_host:'+p] \
                = vector_host_population_contact_sweep_dic[p]

        if len(param_sweep_dic) == 0:
            raise ValueError(
                'param_sweep_dic, host_migration_sweep_dic, vector_migration_sweep_dic, host_host_population_contact_sweep_dic, host_vector_population_contact_sweep_dic, and vector_host_population_contact_sweep_dic cannot all be empty in runParamSweep()'
                )

        params = param_sweep_dic.keys()
        value_lists = [ param_sweep_dic[param] for param in params ]
        combinations = list( it.product( *value_lists ) ) * replicates

        param_df = pd.DataFrame(combinations)
        param_df.columns = params
        results = {}

        print('Starting parallel simulations...')

        def run(param_values):
            model = self.deepCopy()
            for i,param_name in enumerate(params):
                if ':' in param_name:
                    pops = param_name.split(':')[1].split(';')
                    if 'pop_size_host:' in param_name:
                        pop = cp.deepcopy( model.populations[pops[0]] )
                        add_hosts = param_values[i] - len(pop.hosts)
                        if add_hosts < 0:
                            raise ValueError(
                                'Value ' + str(param_values[i]) + ' assigned to ' + pops[0] + ' in host_population_size_sweep must be greater or equal to the population\'s original number of hosts.'
                                )
                        else:
                            pop.addHosts(add_hosts)
                            model.populations[pops[0]] = pop
                            pop.model = model

                    elif 'pop_size_vector:' in param_name:
                        pop = cp.deepcopy( model.populations[pops[0]] )
                        add_vectors = param_values[i] - len(pop.vectors)
                        if add_vectors < 0:
                            raise ValueError(
                                'Values ' + str(param_values[i]) + ' assigned to ' + pops[0] + ' in vector_population_size_sweep must be greater or equal to the population\'s original number of vectors.'
                                )
                        else:
                            pop.addVectors(add_vectors)
                            model.populations[pops[0]] = pop
                            pop.model = model

                    elif 'migrate_host:' in param_name:
                        new_pops = [
                            cp.deepcopy( model.populations[pops[0]] ),
                            cp.deepcopy( model.populations[pops[1]] )
                            ]
                        new_pops[0].model = model
                        new_pops[1].model = model
                        model.populations[pops[0]] = new_pops[0]
                        model.populations[pops[1]] = new_pops[1]
                        model.linkPopulationsHostMigration(
                            new_pops[0],new_pops[1],param_values[i]
                            )
                    elif 'migrate_vector:' in param_name:
                        new_pops = [
                            cp.deepcopy( model.populations[pops[0]] ),
                            cp.deepcopy( model.populations[pops[1]] )
                            ]
                        new_pops[0].model = model
                        new_pops[1].model = model
                        model.populations[pops[0]] = new_pops[0]
                        model.populations[pops[1]] = new_pops[1]
                        model.linkPopulationsVectorMigration(
                            pops[0],pops[1],param_values[i]
                            )
                    elif 'population_contact_host_host:' in param_name:
                        new_pops = [
                            cp.deepcopy( model.populations[pops[0]] ),
                            cp.deepcopy( model.populations[pops[1]] )
                            ]
                        new_pops[0].model = model
                        new_pops[1].model = model
                        model.populations[pops[0]] = new_pops[0]
                        model.populations[pops[1]] = new_pops[1]
                        model.linkPopulationsHostHostContact(
                            pops[0],pops[1],param_values[i]
                            )
                        model.linkPopulationsHostHostContact(
                            pops[1],pops[0],param_values[i]
                            )
                    elif 'population_contact_host_vector:' in param_name:
                        new_pops = [
                            cp.deepcopy( model.populations[pops[0]] ),
                            cp.deepcopy( model.populations[pops[1]] )
                            ]
                        new_pops[0].model = model
                        new_pops[1].model = model
                        model.populations[pops[0]] = new_pops[0]
                        model.populations[pops[1]] = new_pops[1]
                        model.linkPopulationsHostVectorContact(
                            pops[0],pops[1],param_values[i]
                            )
                        model.linkPopulationsHostVectorContact(
                            pops[1],pops[0],param_values[i]
                            )
                    elif 'population_contact_vector_host:' in param_name:
                        new_pops = [
                            cp.deepcopy( model.populations[pops[0]] ),
                            cp.deepcopy( model.populations[pops[1]] )
                            ]
                        new_pops[0].model = model
                        new_pops[1].model = model
                        model.populations[pops[0]] = new_pops[0]
                        model.populations[pops[1]] = new_pops[1]
                        model.linkPopulationsVectorHostContact(
                            pops[0],pops[1],param_values[i]
                            )
                        model.linkPopulationsVectorHostContact(
                            pops[1],pops[0],param_values[i]
                            )
                else:
                    setattr(model.setups[setup_id],param_name,param_values[i])

            for name,pop in model.populations.items():
                pop.setSetup( model.setups[pop.setup.id] )

            sim = Gillespie(model)
            mod = sim.run(
                t0,tf,time_sampling=-1,
                host_sampling=host_sampling,vector_sampling=vector_sampling
                )[tf]
            mod.history = { tf:mod }
            return mod

        return (
            param_df,
            jl.Parallel(n_jobs=n_cores, verbose=10, **kwargs) (
                jl.delayed( run ) (param_values)
                    for param_values in combinations
                )
            )

    def copyState(self,host_sampling=0,vector_sampling=0):
        """Returns a slimmed-down representation of the current model state.

        Keyword arguments:
            host_sampling (int >= 0): how many hosts to skip before saving one in a snapshot
                of the system state (saves all by default). Defaults to 0.
            vector_sampling (int >= 0): how many vectors to skip before saving one in a
                snapshot of the system state (saves all by default). Defaults to 0.

        Returns:
            Model object with current population host and vector lists.
        """

        copy = Model()

        copy.populations = {
            id: p.copyState(host_sampling,vector_sampling)
            for id,p in self.populations.items()
            }

        return copy

    def deepCopy(self):
        """Returns a full copy of the current model with inner references.

        Returns:
            Copied Model object.
        """

        model = cp.deepcopy(self)
        for intervention in model.interventions:
            intervention.model = model
        for pop in model.populations:
            model.populations[pop].model = model
            for h in model.populations[pop].hosts:
                h.population = model.populations[pop]
            for v in model.populations[pop].vectors:
                v.population = model.populations[pop]

        return model


    ### Output and Plots: ###

    def saveToDataFrame(self,save_to_file,n_cores=0,**kwargs):
        """Save status of model to dataframe, write to file location given.

        Creates a pandas Dataframe in long format with the given model history,
        with one host or vector per simulation time in each row, and columns:
            
        - Time - simulation time of entry
        - Population - ID of this host/vector's population
        - Organism - host/vector
        - ID - ID of host/vector
        - Pathogens - all genomes present in this host/vector separated by ;
        - Protection - all genomes present in this host/vector separated by ;
        - Alive - whether host/vector is alive at this time, True/False

        Arguments:
            save_to_file (String): file path and name to save model data under.

        Keyword arguments:
            n_cores (int >= 0): number of cores to parallelize file export across, if 0, all
                cores available are used. Defaults to 0.
            **kwargs: additional arguents for joblib multiprocessing.

        Returns:
            `pandas dataframe` with model history as described above.
        """

        data = saveToDf(
            self.history,save_to_file,n_cores,**kwargs
            )

        return data

    def getPathogens(self, dat, save_to_file=""):
        """Create Dataframe with counts for all pathogen genomes in data.

        Returns sorted pandas Dataframe with counts for occurrences of all
        pathogen genomes in data passed.

        Arguments:
            data (pandas DataFrame): dataframe with model history as produced by `saveToDf` function.

        Keyword arguments:
            save_to_file (String): file path and name to save model data under, no saving
                occurs if empty string. Defaults to "".

        Returns:
            `pandas dataframe` with Series as described above.
        """

        return getPathogens(dat, save_to_file=save_to_file)

    def getProtections(self, dat, save_to_file=""):
        """Create Dataframe with counts for all protection sequences in data.

        Returns sorted `pandas Dataframe` with counts for occurrences of all
        protection sequences in data passed.

        Arguments:
            data (dataframe): dataframe with model history as produced by `saveToDf` function.

        Keyword arguments:
            save_to_file (String): file path and name to save model data under, no saving
                occurs if empty string. Defaults to "".

        Returns:
            pandas DataFrame with Series as described above.
        """

        return getProtections(dat, save_to_file=save_to_file)

    def populationsPlot(
            self, file_name, data, compartment='Infected',
            hosts=True, vectors=False, num_top_populations=7,
            track_specific_populations=[], save_data_to_file="",
            x_label='Time', y_label='Hosts', figsize=(8, 4), dpi=200,
            palette=CB_PALETTE, stacked=False):
        """Create plot with aggregated totals per population across time.

        Creates a line or stacked line plot with dynamics of a compartment
        across populations in the model, with one line for each population.

        A host or vector is considered part of the recovered compartment
        if it has protection sequences of any kind and is not infected.

        Arguments:
            file_name (String): file path, name, and extension to save plot under.
            data (pandas DataFrame): dataframe with model history as produced by saveToDf function.

        Keyword arguments:
            compartment (String): subset of hosts/vectors to count totals of, can be either
                'Naive','Infected','Recovered', or 'Dead'. Defaults to 'Infected'.
            hosts (Boolean): whether to count hosts. Defaults to True.
            vectors (Boolean) whether to count vectors. Defaults to False.
            num_top_populations (int): how many populations to count separately and
                include as columns, remainder will be counted under column "Other";
                if <0, includes all populations in model. Defaults to 7.
            track_specific_populations (list of Strings): contains IDs of specific populations to
                have as a separate column if not part of the top num_top_populations
                populations. Defaults to [].
            save_data_to_file (String): file path and name to save model plot data under,
                no saving occurs if empty string. Defaults to "".
            x_label (String): X axis title. Defaults to 'Time'.
            y_label (String): Y axis title. Defaults to 'Hosts'.
            legend_title (String): legend title. Defaults to 'Population'.
            legend_values (list of Strings): labels for each trace, if empty list, uses population
                IDs. Defaults to [].
            figsize (array-like of two ints): dimensions of figure. Defaults to (8,4).
            dpi (int): figure resolution. Defaults to 200.
            palette (list of color Strings): color palette to use for traces. Defaults to `CB_PALETTE`.
            stacked (Boolean): whether to draw a regular line plot instead of a stacked one. Defaults to False.

        Returns:
            `axis object` for plot with model population dynamics as described above.
        """

        return populationsPlot(
            file_name, data, compartment=compartment, hosts=hosts,
            vectors=vectors, num_top_populations=num_top_populations,
            track_specific_populations=track_specific_populations,
            save_data_to_file=save_data_to_file,
            x_label=x_label, y_label=y_label, figsize=figsize, dpi=dpi,
            palette=palette, stacked=stacked
            )

    def compartmentPlot(
            self, file_name, data, populations=[], hosts=True, vectors=False,
            save_data_to_file="", x_label='Time', y_label='Hosts',
            figsize=(8, 4), dpi=200, palette=CB_PALETTE, stacked=False):
        """Create plot with number of naive,inf,rec,dead hosts/vectors vs. time.

        Creates a line or stacked line plot with dynamics of all compartments
        (naive, infected, recovered, dead) across selected populations in the
        model, with one line for each compartment.

        A host or vector is considered part of the recovered compartment
        if it has protection sequences of any kind and is not infected.

        Arguments:
            file_name (String): file path, name, and extension to save plot under.
            data (pandas DataFrame): dataframe with model history as produced by `saveToDf` function.

        Keyword arguments:
            populations (list of Strings): IDs of populations to include in analysis; if empty, uses
                all populations in model. Defaults to [].
            hosts (Boolean): whether to count hosts. Defaults to True.
            vectors (Boolean): whether to count vectors. Defaults to False.
            save_data_to_file (String): file path and name to save model data under, no
                saving occurs if empty string. Defaults to "".
            x_label (String): X axis title. Defaults to 'Time'.
            y_label (String): Y axis title. Defaults to 'Hosts'.
            legend_title (String): legend title. Defaults to 'Population'.
            legend_values (list of Strings): labels for each trace, if empty list, uses population
                IDs. Defaults to [].
            figsize (array-like of two ints): dimensions of figure. Defaults to (8,4).
            dpi (int)): figure resolution. Defaults to 200.
            palette (list of color Strings): color palette to use for traces. Defaults to `CB_PALETTE`.
            stacked -- whether to draw a regular line plot instead of a stacked one
                (default False, Boolean)

        Returns:
            axis object for plot with model compartment dynamics as described above
        """

        return compartmentPlot(
            file_name, data, populations=populations, hosts=hosts,
            vectors=vectors, save_data_to_file=save_data_to_file,
            x_label=x_label, y_label=y_label, figsize=figsize, dpi=dpi,
            palette=palette, stacked=stacked
            )

    def compositionPlot(
            self, file_name, data, composition_dataframe=None, populations=[],
            type_of_composition='Pathogens', hosts=True, vectors=False,
            num_top_sequences=7, track_specific_sequences=[],
            save_data_to_file="", x_label='Time', y_label='Infections',
            figsize=(8, 4), dpi=200, palette=CB_PALETTE, stacked=True,
            remove_legend=False, genomic_positions=[],population_fraction=False,
            count_individuals_based_on_model=None,
            legend_title='Genotype', legend_values=[], **kwargs):
        """Create plot with counts for pathogen genomes or resistance vs. time.

        Creates a line or stacked line plot with dynamics of the pathogen
        strains or protection sequences across selected populations in the
        model, with one line for each pathogen genome or protection sequence
        being shown.

        Of note: sum of totals for all sequences in one time point does not
        necessarily equal the number of infected hosts and/or vectors, given
        multiple infections in the same host/vector are counted separately.

        Arguments:
            file_name (String): file path, name, and extension to save plot under.
            data (pandas DataFrame): dataframe with model history as produced by `saveToDf` function.

        Keyword arguments:
            composition_dataframe (pandas DataFrame): output of compositionDf() if already computed
                Defaults to None.
            populations (list of Strings): IDs of populations to include in analysis; if empty, uses
                all populations in model. Defaults to [].
            type_of_composition (String) field of data to count totals of, can be either
                'Pathogens' or 'Protection'. Defaults to 'Pathogens'.
            hosts (Boolean) whether to count hosts. Defaults to True.
            vectors (Boolean): whether to count vectors. Defaults to False.
            num_top_sequences (int): how many sequences to count separately and include
                as columns, remainder will be counted under column "Other"; if <0,
                includes all genomes in model. Defaults to 7.
            track_specific_sequences (list of Strings): contains specific sequences to have
                as a separate column if not part of the top num_top_sequences
                sequences. Defaults to [].
            genomic_positions (list of lists of int): list in which each element is a list with loci
                positions to extract (e.g. genomic_positions=[ [0,3], [5,6] ]
                extracts positions 0, 1, 2, and 5 from each genome); if empty, takes
                full genomes. Defaults to [].
            count_individuals_based_on_model (None or Model): `Model` object with populations and
                fitness functions used to evaluate the most fit pathogen genome in
                each host/vector in order to count only a single pathogen per
                host/vector, as opposed to all pathogens within each host/vector; if
                None, counts all pathogens. Defaults to None.
            save_data_to_file (String): file path and name to save model data under, no
                saving occurs if empty string. Defaults to "".
            x_label (String): X axis title. Defaults to 'Time'.
            y_label (String): Y axis title. Defaults to 'Hosts'.
            legend_title (String): legend title. Defaults to 'Population'.
            legend_values (list of Strings): labels for each trace, if empty list, uses population
                IDs. Defaults to [].
            figsize (array-like of two ints): dimensions of figure. Defaults to (8,4).
            dpi (int): figure resolution. Defaults to 200.
            palette (list of color Strings): color palette to use for traces. Defaults to `CB_PALETTE`.
            stacked (Boolean): whether to draw a regular line plot instead of a stacked one. Defaults to False.
            remove_legend (Boolean): whether to print the sequences on the figure legend
                instead of printing them on a separate csv file. Defaults to True.
            **kwargs: additional arguents for joblib multiprocessing.

        Returns:
            axis object for plot with model sequence composition dynamics as described.
        """

        return compositionPlot(
            file_name, data,
            composition_dataframe=composition_dataframe,populations=populations,
            type_of_composition=type_of_composition, hosts=hosts,
            vectors=vectors, num_top_sequences=num_top_sequences,
            track_specific_sequences=track_specific_sequences,
            save_data_to_file=save_data_to_file,
            x_label=x_label, y_label=y_label, figsize=figsize, dpi=dpi,
            palette=palette, stacked=stacked, remove_legend=remove_legend,
            genomic_positions=genomic_positions,
            count_individuals_based_on_model=count_individuals_based_on_model,
            population_fraction=population_fraction, legend_title=legend_title,
            legend_values=legend_values,
            **kwargs
            )

    def clustermap(
            self,
            file_name, data, num_top_sequences=-1, track_specific_sequences=[],
            seq_names=[], n_cores=0, method='weighted', metric='euclidean',
            save_data_to_file="", legend_title='Distance', legend_values=[],
            figsize=(10,10), dpi=200, color_map=DEF_CMAP):
        """Create a heatmap and dendrogram for pathogen genomes in data passed.

        Arguments:
            file_name (String): file path, name, and extension to save plot under.
            data (pandas DataFrame): dataframe with model history as produced by `saveToDf` 
                function.

        Keyword arguments:
            num_top_sequences (int): how many sequences to include in matrix; if <0,
                includes all genomes in data passed. Defaults to -1.
            track_specific_sequences (list of Strings): contains specific sequences to include 
                in matrix if not part of the top num_top_sequences sequences. Defaults to [].
            seq_names (list ofStrings): list with names to be used for sequence labels in matrix
                must be of same length as number of sequences to be displayed; if empty, uses sequences 
                themselves. Defaults to [].
            n_cores (int >= 0): number of cores to parallelize distance compute across, if 0,
                all cores available are used. Defaults to 0.
            method (String): clustering algorithm to use with seaborn clustermap. Defaults to 'weighted'.
            metric (String): distance metric to use with seaborn clustermap. Defaults to 'euclidean'.
            save_data_to_file (String): file path and name to save model data under, no
                saving occurs if empty string. Defaults to "".
            legend_title (String): legend title. Defaults to 'Distance'.
            figsize (array-like of two ints): dimensions of figure. Defaults to (8,4).
            dpi (int): figure resolution. Defaults to 200.
            color_map (matplotlib cmap object): color map to use for traces. Defaults to `DEF_CMAP`.

        Returns:
            figure object for plot with heatmap and dendrogram as described.
        """

        return clustermap(
                file_name, data, num_top_sequences=num_top_sequences,
                track_specific_sequences=track_specific_sequences,
                seq_names=seq_names, n_cores=n_cores, method=method,
                metric=metric, save_data_to_file=save_data_to_file,
                legend_title=legend_title, legend_values=legend_values,
                figsize=figsize, dpi=dpi, color_map=color_map
                )

    def pathogenDistanceHistory(
        self,
        data, samples=-1, num_top_sequences=-1, track_specific_sequences=[],
        seq_names=[], n_cores=0, save_to_file=''):
        """Create DataFrame with pairwise Hamming distances for pathogen
        sequences in data.

        DataFrame has indexes and columns named according to genomes or argument
        seq_names, if passed. Distance is measured as percent Hamming distance
        from an optimal genome sequence.

        Arguments:
            data (pandas DataFrame): dataframe with model history as produced by `saveToDf` 
                function.

        Keyword arguments:
            samples (int): how many timepoints to uniformly sample from the total
                timecourse; if <0, takes all timepoints. Defaults to 1.
            num_top_sequences (int) how many sequences to include in matrix; if <0,
                includes all genomes in data passed. Defaults to -1.
            track_specific_sequences (list of Strings): contains specific sequences to include in
                matrix if not part of the top num_top_sequences sequences. Defaults to [].
            seq_names (list of Strings): list with names to be used for sequence labels in matrix
                must be of same length as number of sequences to be displayed; if
                empty, uses sequences themselves. Defaults to [].
            n_cores (int >= 0): number of cores to parallelize distance compute across, if 0,
                all cores available are used. Defaults to 0.
            save_to_file (String): file path and name to save model data under, no saving
                occurs if empty string. Defaults to "".

        Returns:
            pandas DataFrame with distance matrix as described above.
        """
        return getPathogenDistanceHistoryDf(data,
            samples=samples, num_top_sequences=num_top_sequences,
            track_specific_sequences=track_specific_sequences,
            seq_names=seq_names, n_cores=n_cores, save_to_file=save_to_file)

    def getGenomeTimes(
        self,
        data, samples=-1, num_top_sequences=-1, track_specific_sequences=[],
        seq_names=[], n_cores=0, save_to_file=''):
        """Create DataFrame with times genomes first appeared during simulation.

        Arguments:
            data (pandas DataFrame): dataframe with model history as produced by `saveToDf` 
                function.

        Keyword arguments:
            samples (int): how many timepoints to uniformly sample from the total
                timecourse; if <0, takes all timepoints. Defaults to 1.
            save_to_file (String): file path and name to save model data under, no saving
                occurs if empty string. Defaults to "".
            n_cores (int): number of cores to parallelize across, if 0, all cores
                available are used. Defaults to 0.

        Returns:
            pandas DataFrame with genomes and times as described above.
        """
        return getGenomeTimesDf(data,
            samples=samples, num_top_sequences=num_top_sequences,
            track_specific_sequences=track_specific_sequences,
            seq_names=seq_names, n_cores=n_cores, save_to_file=save_to_file)


    def getCompositionData(
            self, data=None, populations=[], type_of_composition='Pathogens',
            hosts=True, vectors=False, num_top_sequences=-1,
            track_specific_sequences=[], genomic_positions=[],
            count_individuals_based_on_model=None, save_data_to_file="",
            n_cores=0, **kwargs):
        """Create dataframe with counts for pathogen genomes or resistance.

        Creates a pandas Dataframe with dynamics of the pathogen strains or
        protection sequences across selected populations in the model,
        with one time point in each row and columns for pathogen genomes or
        protection sequences.

        Of note: sum of totals for all sequences in one time point does not
        necessarily equal the number of infected hosts and/or vectors, given
        multiple infections in the same host/vector are counted separately.

        Keyword arguments:
            data (pandas DataFrame): dataframe with model history as produced by `saveToDf` 
                function; if None, computes this dataframe and saves it under 'raw_data_'+'save_data_to_file'. 
                Defaults to None.
            populations (list of Strings): IDs of populations to include in analysis; if empty, uses
                all populations in model. Defaults to [].
            type_of_composition (String): field of data to count totals of, can be either
                'Pathogens' or 'Protection'. Defaults to 'Pathogens'.
            hosts (Boolean): whether to count hosts. Defaults to True.
            vectors (Boolean): whether to count vectors. Defaults to False.
            num_top_sequences (int): how many sequences to count separately and include
                as columns, remainder will be counted under column "Other"; if <0,
                includes all genomes in model. Defaults to -1.
            track_specific_sequences (list of Strings): contains specific sequences to have
                as a separate column if not part of the top num_top_sequences
                sequences. Defaults to [].
            genomic_positions (list of lists of int): list in which each element is a list with 
                loci positions to extract (e.g. genomic_positions=[ [0,3], [5,6] ]
                extracts positions 0, 1, 2, and 5 from each genome); if empty, takes
                full genomes. Defaults to [].
            count_individuals_based_on_model (None or Model object): Model object with populations and
                fitness functions used to evaluate the most fit pathogen genome in
                each host/vector in order to count only a single pathogen per
                host/vector, asopposed to all pathogens within each host/vector; if
                None, counts all pathogens. Defaults to None.
            save_data_to_file (String): file path and name to save model data under, no
                saving occurs if empty string. Defaults to "".
            n_cores (int): number of cores to parallelize processing across, if 0, all
                cores available are used. Defaults to 0.
            **kwargs: additional arguents for joblib multiprocessing.

        Returns:
            pandas DataFrame with model sequence composition dynamics as described
                above.
        """

        if data is None:
            data = saveToDf(
                self.history,'raw_data_'+save_to_file,n_cores,verbose=verbose,
                **kwargs
                )

        return compositionDf(
            data, populations=populations,
            type_of_composition=type_of_composition,
            hosts=hosts, vectors=vectors, num_top_sequences=num_top_sequences,
            track_specific_sequences=track_specific_sequences,
            genomic_positions=genomic_positions,
            count_individuals_based_on_model=count_individuals_based_on_model,
            save_to_file=save_data_to_file, n_cores=n_cores, **kwargs
            )

    ### Model interventions: ###

    def newPopulation(self, id, setup_name, num_hosts=0, num_vectors=0):
        """Create a new Population object with setup parameters.

        If population ID is already in use, appends _2 to it

        Arguments:
            id (String): unique identifier for this population in the model.
            setup_name (Setup object): setup object with parameters for this population.

        Keyword arguments:
            num_hosts (int >= 0): number of hosts to initialize population with. Defaults to 100.
            num_vectors (int >= 0): number of vectors to initialize population with. Defaults to 100.
        """

        if id in self.populations.keys():
            id = id+'_2'

        self.populations[id] = Population(
            self, id, self.setups[setup_name], num_hosts, num_vectors
            )

        for p in self.populations:
            self.populations[id].setHostMigrationNeighbor(self.populations[p],0)
            self.populations[id].setVectorMigrationNeighbor(
                self.populations[p],0
                )
            self.populations[id].setHostHostPopulationContactNeighbor(
                self.populations[p],0
                )
            self.populations[id].setVectorHostPopulationContactNeighbor(
                self.populations[p],0
                )
            self.populations[id].setHostVectorPopulationContactNeighbor(
                self.populations[p],0
                )

            self.populations[p].setHostMigrationNeighbor(
                self.populations[id],0
                )
            self.populations[p].setVectorMigrationNeighbor(
                self.populations[id],0
                )
            self.populations[p].setHostHostPopulationContactNeighbor(
                self.populations[id],0
                )
            self.populations[p].setVectorHostPopulationContactNeighbor(
                self.populations[id],0
                )
            self.populations[p].setHostVectorPopulationContactNeighbor(
                self.populations[id],0
                )

    def linkPopulationsHostMigration(self, pop1_id, pop2_id, rate):
        """Set host migration rate from one population towards another.

        Arguments:
            pop1_id (String): origin population for which migration rate will be specified.
            pop1_id (String): destination population for which migration rate will be
                specified.
            rate (number >= 0): migration rate from one population to the neighbor; evts/time.
        """

        self.populations[pop1_id].setHostMigrationNeighbor(
            self.populations[pop2_id], rate
            )

    def linkPopulationsVectorMigration(self, pop1_id, pop2_id, rate):
        """Set vector migration rate from one population towards another.

        Arguments:
            pop1_id (String): origin population for which migration rate will be specified.
            pop1_id (String): destination population for which migration rate will be
                specified.
            rate (number >= 0): migration rate from one population to the neighbor; evts/time.
        """

        self.populations[pop1_id].setVectorMigrationNeighbor(
            self.populations[pop2_id], rate
            )

    def linkPopulationsHostHostContact(self, pop1_id, pop2_id, rate):
        """Set host-host contact rate from one population towards another.

        Arguments:
            pop1_id (String): origin population for which inter-population contact rate
                will be specified.
            pop1_id (String): destination population for which inter-population contact
                rate will be specified.
            rate (number >= 0): inter-population contact rate from one population to the
                neighbor; evts/time.
        """

        self.populations[pop1_id].setHostHostPopulationContactNeighbor(
            self.populations[pop2_id], rate
            )

    def linkPopulationsHostVectorContact(self, pop1_id, pop2_id, rate):
        """Set host-vector contact rate from one population towards another.

        Arguments:
            pop1_id (String): origin population for which inter-population contact rate
                will be specified.
            pop1_id (String): destination population for which inter-population contact
                rate will be specified.
            rate (number >= 0): inter-population contact rate from one population to the
                neighbor; evts/time.
        """

        self.populations[pop1_id].setHostVectorPopulationContactNeighbor(
            self.populations[pop2_id], rate
            )

    def linkPopulationsVectorHostContact(self, pop1_id, pop2_id, rate):
        """Set vector-host contact rate from one population towards another.

        Arguments:
            pop1_id (String): origin population for which inter-population contact rate
                will be specified.
            pop1_id (String): destination population for which inter-population contact
                rate will be specified.
            rate (number >= 0): inter-population contact rate from one population to the
                neighbor; evts/time.
        """

        self.populations[pop1_id].setVectorHostPopulationContactNeighbor(
            self.populations[pop2_id], rate
            )

    def createInterconnectedPopulations(
            self, num_populations, id_prefix, setup_name,
            host_migration_rate=0, vector_migration_rate=0,
            host_host_contact_rate=0,
            host_vector_contact_rate=0, vector_host_contact_rate=0,
            num_hosts=100, num_vectors=100):
        """Create new populations, link all of them to each other.

        All populations in this cluster are linked with the same migration rate,
        starting number of hosts and vectors, and setup parameters. Their IDs
        are numbered onto prefix given as 'id_prefix_0', 'id_prefix_1',
        'id_prefix_2', etc.

        Arguments:
            num_populations (int): number of populations to be created.
            id_prefix (String): prefix for IDs to be used for this population in the model.
            setup_name (Setup object): setup object with parameters for all populations.

        Keyword arguments:
            host_migration_rate (number >= 0): host migration rate between populations;
                evts/time. Defaults to 0.
            vector_migration_rate (number >= 0): vector migration rate between populations;
                evts/time. Defaults to 0.
            host_host_contact_rate (number >= 0): host-host inter-population contact rate
                between populations; evts/time. Defaults to 0.
            host_vector_contact_rate (number >= 0): host-vector inter-population contact rate
                between populations; evts/time. Defaults to 0.
            vector_host_contact_rate (number >= 0): vector-host inter-population contact rate
                between populations; evts/time. Defaults to 0.
            num_hosts (int): number of hosts to initialize population with. Defaults to 100.
            num_vectors (int): number of hosts to initialize population with. Defaults to 100.
        """

        new_pops = [
            Population(
                self, str(id_prefix) + str(i), self.setups[setup_name],
                num_hosts, num_vectors
                ) for i in range(num_populations)
            ]
        new_pop_ids = []
        for pop in new_pops:
            if pop.id in self.populations.keys():
                pop.id = pop.id+'_2'

            self.populations[pop.id] = pop
            new_pop_ids.append(pop.id)

            for p in self.populations:
                pop.setHostMigrationNeighbor(self.populations[p],0)
                pop.setVectorMigrationNeighbor(self.populations[p],0)
                pop.setHostHostPopulationContactNeighbor(self.populations[p],0)
                pop.setHostVectorPopulationContactNeighbor(self.populations[p],0)
                pop.setVectorHostPopulationContactNeighbor(self.populations[p],0)

                self.populations[p].setHostMigrationNeighbor(pop,0)
                self.populations[p].setVectorMigrationNeighbor(pop,0)
                self.populations[p].setHostHostPopulationContactNeighbor(pop,0)
                self.populations[p].setHostVectorPopulationContactNeighbor(pop,0)
                self.populations[p].setVectorHostPopulationContactNeighbor(pop,0)

        for p1_id in new_pop_ids:
            for p2_id in new_pop_ids:
                self.linkPopulationsHostMigration(
                    p1_id,p2_id,host_migration_rate
                    )
                self.linkPopulationsVectorMigration(
                    p1_id,p2_id,vector_migration_rate
                    )
                self.linkPopulationsHostHostContact(
                    p1_id,p2_id,host_host_contact_rate
                    )
                self.linkPopulationsHostVectorContact(
                    p1_id,p2_id,host_vector_contact_rate
                    )
                self.linkPopulationsVectorHostContact(
                    p1_id,p2_id,vector_host_contact_rate
                    )

    def newHostGroup(self, pop_id, group_id, hosts=-1, type='any'):
        """Return a list of random hosts in population.

        Arguments:
            pop_id (String): ID of population to be sampled from.
            group_id (String): ID to name group with.

        Keyword arguments:
            hosts (number): number of hosts to be sampled randomly: if <0, samples from
                whole population; if <1, takes that fraction of population; if >=1,
                samples that integer number of hosts. Defaults to -1.
            type (String = {'healthy', 'infected', 'any'}): whether to sample healthy hosts only, 
                infected hosts only, or any hosts. Defaults to 'any'.

        Returns:
            list containing sampled hosts.
        """

        self.groups[group_id] = self.populations[pop_id].newHostGroup(
            hosts, type
            )

    def newVectorGroup(self, pop_id, group_id, vectors=-1, type='any'):
        """Return a list of random vectors in population.

        Arguments:
            pop_id (String): ID of population to be sampled from.
            group_id (String): ID to name group with.

        Keyword arguments:
            vectors (number): number of vectors to be sampled randomly: if <0, samples from
                whole population; if <1, takes that fraction of population; if >=1,
                samples that integer number of vectors. Defaults to -1.
            type (String = {'healthy', 'infected', 'any'}): whether to sample healthy vectors only, infected vectors
                only, or any vectors. Defaults to 'any'.

        Returns:
            list containing sampled vectors.
        """

        self.groups[group_id] = self.populations[pop_id].newVectorGroup(
            vectors, type
            )

    def addHosts(self, pop_id, num_hosts):
        """Add a number of healthy hosts to population, return list with them.

        Arguments:
            pop_id (String): ID of population to be modified.
            num_hosts (int): number of hosts to be added.

        Returns:
            list containing new hosts.
        """

        self.populations[pop_id].addHosts(num_hosts)

    def addVectors(self, pop_id, num_vectors):
        """Add a number of healthy vectors to population, return list with them.

        Arguments:
            pop_id (String): ID of population to be modified.
            num_vectors (int): number of vectors to be added.

        Returns:
            list containing new vectors.
        """

        self.populations[pop_id].addVectors(num_vectors)

    def removeHosts(self, pop_id, num_hosts_or_list):
        """Remove a number of specified or random hosts from population.

        Arguments:
            pop_id (String): ID of population to be modified.
            num_hosts_or_list (int or list of Hosts): number of hosts to be sampled randomly for removal
                or list of hosts to be removed, must be hosts in this population.
        """

        self.populations[pop_id].removeHosts(num_hosts_or_list)

    def removeVectors(self, pop_id, num_vectors_or_list):
        """Remove a number of specified or random vectors from population.

        Arguments:
            pop_id (String): ID of population to be modified.
            num_vectors_or_list (int or list of Vectors): number of vectors to be sampled randomly for
                removal or list of vectors to be removed, must be vectors in this
                population.
        """

        self.populations[pop_id].removeVectors(num_vectors_or_list)

    def addPathogensToHosts(self, pop_id, genomes_numbers, group_id=""):
        """Add specified pathogens to random hosts, optionally from a list.

        Arguments:
            pop_id (String): ID of population to be modified.
            genomes_numbers (dict with keys=Strings, values=int) dictionary containing pathogen 
                genomes to add as keys and number of hosts each one will be added to as values.

        Keyword arguments:
            group_id (String): ID of specific hosts to sample from, if empty, samples
                from whole population. Defaults to "".
        """

        if group_id == "":
            hosts = self.populations[pop_id].hosts
        else:
            hosts = self.groups[group_id]

        self.populations[pop_id].addPathogensToHosts(genomes_numbers,hosts)

    def addPathogensToVectors(self, pop_id, genomes_numbers, group_id=""):
        """Add specified pathogens to random vectors, optionally from a list.

        Arguments:
            pop_id (String): ID of population to be modified.
            genomes_numbers (dict with keys=Strings, values=int): dictionary containing pathogen 
                genomes to add as keys and number of vectors each one will be added to as values.

        Keyword arguments:
            group_id (String): ID of specific vectors to sample from, if empty, samples
                from whole population. Defaults to "".
        """

        if group_id == "":
            vectors = self.populations[pop_id].vectors
        else:
            vectors = self.groups[group_id]

        self.populations[pop_id].addPathogensToVectors(genomes_numbers,vectors)

    def treatHosts(self, pop_id, frac_hosts, resistance_seqs, group_id=""):
        """Treat random fraction of infected hosts against some infection.

        Removes all infections with genotypes susceptible to given treatment.
        Pathogens are removed if they are missing at least one of the sequences
        in resistance_seqs from their genome. Removes this organism from
        population infected list and adds to healthy list if appropriate.

        Arguments:
            pop_id (String): ID of population to be modified.
            frac_hosts (number between 0 and 1): fraction of hosts considered to be randomly selected.
            resistance_seqs (list of Strings): contains sequences required for treatment resistance.

        Keyword arguments:
            group_id (String): ID of specific hosts to sample from, if empty, samples
                from whole population. Defaults to "".
        """

        if group_id == "":
            hosts = self.populations[pop_id].hosts
        else:
            hosts = self.groups[group_id]

        self.populations[pop_id].treatHosts(frac_hosts,resistance_seqs,hosts)

    def treatVectors(self, pop_id, frac_vectors, resistance_seqs, group_id=""):
        """Treat random fraction of infected vectors agains some infection.

        Removes all infections with genotypes susceptible to given treatment.
        Pathogens are removed if they are missing at least one of the sequences
        in resistance_seqs from their genome. Removes this organism from
        population infected list and adds to healthy list if appropriate.

        Arguments:
            pop_id (String): ID of population to be modified.
            frac_vectors (number between 0 and 1): fraction of vectors considered to be 
                randomly selected.
            resistance_seqs (list of Strings): contains sequences required for treatment resistance.

        Keyword arguments:
            group_id (String): ID of specific vectors to sample from, if empty, samples
                from whole population. Defaults to "".
        """

        if group_id == "":
            vectors = self.populations[pop_id].vectors
        else:
            vectors = self.groups[group_id]

        self.populations[pop_id].treatVectors(
            frac_vectors,resistance_seqs,vectors
            )

    def protectHosts(
            self, pop_id, frac_hosts, protection_sequence, group_id=""):
        """Protect a random fraction of infected hosts against some infection.

        Adds protection sequence specified to a random fraction of the hosts
        specified. Does not cure them if they are already infected.

        Arguments:
            pop_id (String): ID of population to be modified.
            frac_hosts (number between 0 and 1): fraction of hosts considered to be 
                randomly selected.
            protection_sequence (String): sequence against which to protect.

        Keyword arguments:
            group_id (String): ID of specific hosts to sample from, if empty, samples
                from whole population. Defaults to "".
        """

        if group_id == "":
            hosts = self.populations[pop_id].hosts
        else:
            hosts = self.groups[group_id]

        self.populations[pop_id].protectHosts(
            frac_hosts,protection_sequence,hosts
            )

    def protectVectors(
            self, pop_id, frac_vectors, protection_sequence, group_id=""):
        """Protect a random fraction of infected vectors against some infection.

        Adds protection sequence specified to a random fraction of the vectors
        specified. Does not cure them if they are already infected.

        Arguments:
            pop_id (String): ID of population to be modified.
            frac_vectors (number between 0 and 1): fraction of vectors considered to be randomly selected.
            protection_sequence (String): sequence against which to protect.

        Keyword arguments:
            group_id (String): ID of specific vectors to sample from, if empty, samples
                from whole population. Defaults to "".
        """

        if group_id == "":
            vectors = self.populations[pop_id].vectors
        else:
            vectors = self.groups[group_id]

        self.populations[pop_id].protectVectors(
            frac_vectors,protection_sequence,vectors
            )

    def wipeProtectionHosts(self, pop_id, group_id=""):
        """Removes all protection sequences from hosts.

        Arguments:
            pop_id (String): ID of population to be modified.

        Keyword arguments:
            group_id (String): ID of specific hosts to sample from, if empty, samples from
                whole from whole population. Defaults to "".
        """

        if group_id == "":
            hosts = self.populations[pop_id].hosts
        else:
            hosts = self.groups[group_id]

        self.populations[pop_id].wipeProtectionHosts(hosts)

    def wipeProtectionVectors(self, pop_id, group_id=""):
        """Removes all protection sequences from vectors.

        Arguments:
            pop_id (String): ID of population to be modified.

        Keyword arguments:
            group_id (String): ID of specific vectors to sample from, if empty, samples
                from whole population. Defaults to "".
        """

        if group_id == "":
            vectors = self.populations[pop_id].vectors
        else:
            vectors = self.groups[group_id]

        self.populations[pop_id].wipeProtectionVectors(vectors)

    ### Modify population parameters: ###

    def setSetup(self, pop_id, setup_id):
        """Assign parameters stored in Setup object to this population.

        Arguments:
            pop_id (String): ID of population to be modified.
            setup_id (String): ID of setup to be assigned.
        """

        self.populations[pop_id].setSetup( self.setups[setup_id] )

    ### Utility: ###

    def customModelFunction(self, function):
        """Returns output of given function, passing this model as a parameter.

        Arguments:
            function (callable): function to be evaluated; must take a Model object as the
                        only parameter.

        Returns:
            output of function passed as parameter.
        """

        return function(self)

    ### Preset fitness functions: ###

    @staticmethod
    def peakLandscape(genome, peak_genome, min_value):
        """Return genome phenotype by decreasing with distance from optimal seq.

        A purifying selection fitness function based on exponential decay of
        fitness as genomes move away from the optimal sequence. Distance is
        measured as percent Hamming distance from an optimal genome sequence.

        Arguments:
            genome (String): the genome to be evaluated.
            peak_genome (String): the genome sequence to measure distance against, has
                value of 1.
            min_value (number 0-1): minimum value at maximum distance from optimal
                genome.

        Return:
            value of genome (number).
        """

        distance = td.hamming(genome, peak_genome) / len(genome)
        value = np.exp( np.log( min_value ) * distance )

        return value

    @staticmethod
    def valleyLandscape(genome, valley_genome, min_value):
        """Return genome phenotype by increasing with distance from worst seq.

        A disruptive selection fitness function based on exponential decay of
        fitness as genomes move closer to the worst possible sequence. Distance
        is measured as percent Hamming distance from the worst possible genome
        sequence.

        Arguments:
            genome (String): the genome to be evaluated.
            valley_genome (String): the genome sequence to measure distance against, has
                value of min_value.
            min_value (number 0-1): fitness value of worst possible genome.

        Return:
            value of genome (number).
        """

        distance = td.hamming(genome, valley_genome) / len(genome)
        value = np.exp( np.log( min_value ) * ( 1 - distance ) )

        return value