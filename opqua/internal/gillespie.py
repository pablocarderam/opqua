
"""Contains class Population."""

import copy as cp
import pandas as pd
import numpy as np
import joblib as jl

from opqua.model import *

# Event ID constants:

MIGRATE = 0
CONTACT_INFECTED_HOST_ANY_HOST = 1
CONTACT_INFECTED_HOST_ANY_VECTOR = 2
CONTACT_HEALTHY_HOST_INFECTED_VECTOR = 3
RECOVER_HOST = 4
RECOVER_VECTOR = 5
MUTATE_HOST = 6
MUTATE_VECTOR = 7
RECOMBINE_HOST = 8
RECOMBINE_VECTOR = 9
KILL_HOST = 10
KILL_VECTOR = 11

class Gillespie(object):
    """Class contains methods for simulating model with Gillespie algorithm.

    Class defines a model's events and methods for changing system state
    according to the possible events and simulating a timecourse using the
    Gillespie algorithm.

    Methods:
    getRates -- returns array containing rates for each event for a given system
        state
    doAction -- carries out an event, modifying system state
    run -- simulates model for a specified length of time
    """

    def __init__(self, model):
        """Create a new Gillespie simulation object.

        Arguments:
        model -- the model this simulation belongs to (Model)
        """

        super(Gillespie, self).__init__() # initialize as parent class object

        # Event IDs
        self.evt_IDs = [
            MIGRATE, CONTACT_INFECTED_HOST_ANY_HOST,
            CONTACT_INFECTED_HOST_ANY_VECTOR,
            CONTACT_HEALTHY_HOST_INFECTED_VECTOR,
            RECOVER_HOST, RECOVER_VECTOR,
            MUTATE_HOST, MUTATE_VECTOR,
            RECOMBINE_HOST, RECOMBINE_VECTOR,
            KILL_HOST, KILL_VECTOR
            ]
            # event IDs in specific order to be used

        self.model = model

    def getRates(self,population_ids):
        """Wrapper for calculating event rates as per current system state.

        Arguments:
        population_ids -- list with ids for every population in the model
            (list of Strings)

        Returns:
        Matrix with rates as values for events (rows) and populations (columns).
        Populations in order given in argument.
        """

        hosts,infected_hosts,healthy_hosts,vectors,infected_vectors, \
        total_migration_rate, \
        contact_rate_host_host,contact_rate_host_vector, \
        recovery_rate_host,recovery_rate_vector, \
        mutate_in_host,mutate_in_vector, \
        recombine_in_host,recombine_in_vector, \
        death_rate_host,death_rate_vector = (
            [0] * len(population_ids), [0] * len(population_ids),
            [0] * len(population_ids), [0] * len(population_ids),
            [0] * len(population_ids), [0] * len(population_ids),
            [0] * len(population_ids), [0] * len(population_ids),
            [0] * len(population_ids), [0] * len(population_ids),
            [0] * len(population_ids), [0] * len(population_ids),
            [0] * len(population_ids), [0] * len(population_ids),
            [0] * len(population_ids), [0] * len(population_ids)
            )

        for i,id in enumerate(population_ids):
            hosts[i] = len( self.model.populations[id].hosts )
            infected_hosts[i] = len( self.model.populations[id].infected_hosts )
            healthy_hosts[i] = len( self.model.populations[id].healthy_hosts )
            vectors[i] = len( self.model.populations[id].vectors )
            infected_vectors[i] = len( self.model.populations[id].infected_vectors )
            total_migration_rate[i] = self.model.populations[id].total_migration_rate
            contact_rate_host_host[i] = self.model.populations[id].contact_rate_host_host
            contact_rate_host_vector[i] = self.model.populations[id].contact_rate_host_vector
            recovery_rate_host[i] = self.model.populations[id].recovery_rate_host
            recovery_rate_vector[i] = self.model.populations[id].recovery_rate_vector
            mutate_in_host[i] = self.model.populations[id].mutate_in_host
            mutate_in_vector[i] = self.model.populations[id].mutate_in_vector
            recombine_in_host[i] = self.model.populations[id].recombine_in_host
            recombine_in_vector[i] = self.model.populations[id].recombine_in_vector
            death_rate_host[i] = self.model.populations[id].death_rate_host
            death_rate_vector[i] = self.model.populations[id].death_rate_vector

        return getRatesExternal( len(self.evt_IDs), len(population_ids),
            np.array(hosts),np.array(infected_hosts),np.array(healthy_hosts),np.array(vectors),np.array(infected_vectors),
            np.array(total_migration_rate),
            np.array(contact_rate_host_host),np.array(contact_rate_host_vector),
            np.array(recovery_rate_host),np.array(recovery_rate_vector),
            np.array(mutate_in_host),np.array(mutate_in_vector),
            np.array(recombine_in_host),np.array(recombine_in_vector),
            np.array(death_rate_host),np.array(death_rate_vector)
            )

    def doAction(self,act,pop,rand):

        """Change system state according to act argument passed

        Arguments:
        act -- defines action to be taken, one of the event ID constants (int)
        pop -- population action will happen in (Population)
        rand -- random number used to define event (number 0-1)

        Returns:
        whether or not the model has changed state (Boolean)
        """

        changed = False

        if act == MIGRATE:
            rand = rand * pop.total_migration_rate
            r_cum = 0
            for neighbor in pop.neighbors:
                r_cum += pop.neighbors[neighbor]
                if r_cum > rand:
                    pop.migrate(neighbor,1,0)
                    changed = True
                    break

        elif act == CONTACT_INFECTED_HOST_ANY_HOST:
            rand = rand * len(pop.infected_hosts)
            infected_host = int( np.floor(rand) )
            other_host = int(
                np.floor( ( rand - infected_host ) * len(pop.hosts) )
                )
            changed = pop.contactInfectedHostAnyHost(infected_host,other_host)

        elif act == CONTACT_INFECTED_HOST_ANY_VECTOR:
            rand = rand * len(pop.infected_hosts)
            infected_host = int( np.floor(rand) )
            vector = int(
                np.floor( ( rand - infected_host ) * len(pop.vectors) )
                )
            changed = pop.contactInfectedHostAnyVector(infected_host,vector)

        elif act == CONTACT_HEALTHY_HOST_INFECTED_VECTOR:
            rand = rand * len(pop.healthy_hosts)
            healthy_host = int( np.floor(rand) )
            infected_vector = int(
                np.floor( ( rand - healthy_host ) * len(pop.infected_vectors) )
                )
            changed = pop.contactHealthyHostInfectedVector(
                healthy_host,infected_vector
                )

        elif act == RECOVER_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            pop.recoverHost(host)
            changed = True

        elif act == RECOVER_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.recoverVector(vector)
            changed = True

        elif act == MUTATE_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            pop.mutateHost(host)
            changed = True

        elif act == MUTATE_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.mutateVector(vector)
            changed = True

        elif act == RECOMBINE_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            pop.recombineHost(host)
            changed = True

        elif act == RECOMBINE_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.recombineVector(vector)
            changed = True

        elif act == KILL_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            rand = rand * len(pop.infected_hosts) - host
            pop.killHost(host,rand)
            changed = True

        elif act == KILL_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            rand = rand * len(pop.infected_vectors) - vector
            pop.killVector(vector,rand)
            changed = True

        return changed

    def run(self,t0,tf,time_sampling=0,host_sampling=0,vector_sampling=0,
            print_every_n_events=1000):

        """Simulate model for a specified time between two time points.

        Simulates a time series using the Gillespie algorithm.

        Arguments:
        t0 -- initial time point to start simulation at (number)
        tf -- initial time point to end simulation at (number)
        time_sampling -- how many events to skip before saving a snapshot of the
            system state (saves all by default), if <0, saves only final state
            (int, default 0)
        host_sampling -- how many hosts to skip before saving one in a snapshot
            of the system state (saves all by default) (int, default 0)
        vector_sampling -- how many vectors to skip before saving one in a
            snapshot of the system state (saves all by default) (int, default 0)
        print_every_n_events -- number of events a message is printed to console
            (int>0, default 1000)

        Returns:
        dictionary containing model state history, with keys=times and
            values=Model objects with model snapshot at that time point
        """

        # Simulation variables
        t_var = t0 # keeps track of time
        history = { 0: cp.deepcopy(self.model) }
        intervention_tracker = 0
            # keeps track of what the next intervention should be
        self.model.interventions = sorted(
            self.model.interventions, key=lambda i: i.time
            )

        print_counter = 0 # only used to track when to print
        sampling_counter = 0 # used to track when to save a snapshot

        while t_var < tf: # repeat until t reaches end of timecourse
            population_ids = list( self.model.populations.keys() )
            r = self.getRates(population_ids) # get event rates in this state
            r_tot = np.sum(r) # sum of all rates

            # Time handling
            if r_tot > 0:
                dt = np.random.exponential( 1/r_tot ) # time until next event
                t_var += dt # add time step to main timer

                if (intervention_tracker < len(self.model.interventions)
                    and t_var
                    >= self.model.interventions[intervention_tracker].time):
                    # if there are any interventions left and if it is time
                    # to make one,
                    self.model.interventions[
                        intervention_tracker
                        ].doIntervention()
                    t_var = self.model.interventions[intervention_tracker].time
                    intervention_tracker += 1 # advance the tracker

                    population_ids = list( self.model.populations.keys() )
                    r = self.getRates(population_ids)
                        # get event rates in this state
                    r_tot = np.sum(r) # sum of all rates
                    dt = np.random.exponential( 1/r_tot )
                        # time until next event
                    t_var += dt # add time step to main timer

                # Event handling
                if t_var < tf: # if still within max time
                    u = np.random.random() * r_tot
                        # random uniform number between 0 and total rate
                    r_cum = 0 # cumulative rate
                    for e in range(r.shape[0]): # for every possible event,
                        for p in range(r.shape[1]):
                            # for every possible population,
                            r_cum += r[e,p]
                                # add this event's rate to cumulative rate
                            if u < r_cum:
                                # if random number is under cumulative rate

                                # print every n events
                                print_counter += 1
                                if print_counter == print_every_n_events:
                                    print_counter = 0
                                    print(
                                        'Simulating time: '
                                        + str(t_var) + ', event ID: ' + str(e)
                                        )

                                changed = self.doAction(
                                    e, self.model.populations[
                                        population_ids[p]
                                        ], ( u - r_cum + r[e,p] ) / r[e,p]
                                    ) # do corresponding action,
                                      # feed in renormalized random number
                                if changed and time_sampling >= 0:
                                        # if state changed and saving history,
                                        # saves history at correct intervals
                                    sampling_counter += 1
                                    if sampling_counter > time_sampling:
                                        sampling_counter = 0
                                        history[t_var] = self.model.copyState()

                                break # exit event loop

                        else: # if the inner loop wasn't broken,
                            continue # continue outer loop

                        break # otherwise, break outer loop
            else: # if at end of time course,
                if intervention_tracker < len(self.model.interventions):
                    print( 'Simulating time: ' + str(t_var), 'END')
                    self.model.interventions[
                        intervention_tracker
                        ].doIntervention()
                    t_var = self.model.interventions[intervention_tracker].time
                    intervention_tracker += 1 # advance the tracker
                else:
                    print( 'Simulating time: ' + str(t_var), 'END')
                    t_var = tf

        history[tf] = cp.deepcopy(self.model)
        history[tf].history = None

        return history


# @jit(nopython=True) testing out possible Numba implementation
def getRatesExternal(num_evt_IDs,num_population_ids,
        hosts,infected_hosts,healthy_hosts,vectors,infected_vectors,
        total_migration_rate,contact_rate_host_host,contact_rate_host_vector,
        recovery_rate_host,recovery_rate_vector,
        mutate_in_host,mutate_in_vector,recombine_in_host,recombine_in_vector,
        death_rate_host,death_rate_vector
        ):
    """Calculates event rates according to current system state.

    Arguments:
    population_ids -- list with ids for every population in the model
        (list of Strings)

    Returns:
    dictionary with event ID constants as keys and rates as values. Includes
        total rate under 'tot' key.
    """

    rates = np.zeros( [ num_evt_IDs, num_population_ids ] )
        # rate array size of event space

    rates[MIGRATE,:] = np.multiply( hosts, total_migration_rate)

    rates[CONTACT_INFECTED_HOST_ANY_HOST,:] = np.multiply(
        np.divide( infected_hosts, np.maximum( hosts, 1 ) ),
        contact_rate_host_host
        )
        # contact rates assumes scaling area--large populations are equally
        # dense as small ones, so contact is constant with both host and
        # vector populations. If you don't want this to happen, modify the
        # population's contact rate accordingly.
        # Examines contacts between infected hosts and all hosts

    rates[CONTACT_INFECTED_HOST_ANY_VECTOR,:] = np.multiply(
        np.divide( infected_hosts, np.maximum( hosts, 1 ) ),
        contact_rate_host_vector
        )
        # contact rates assumes scaling area--large populations are equally
        # dense as small ones, so contact is constant with both host and
        # vector populations. If you don't want this to happen, modify the
        # population's contact rate accordingly.
        # Examines contacts between infected hosts and all vectors

    rates[CONTACT_HEALTHY_HOST_INFECTED_VECTOR,:] = np.multiply(
        np.multiply(
            np.divide( healthy_hosts, np.maximum( hosts, 1 ) ),
            np.divide( infected_vectors, np.maximum( vectors, 1 ) )
            ),
        contact_rate_host_vector
        )
        # contact rates assumes scaling area--large populations are equally
        # dense as small ones, so contact is constant with both host and
        # vector populations. If you don't want this to happen, modify the
        # population's contact rate accordingly.
        # Examines contacts between healthy hosts and infected vectors
        # (infected-infected contacts are considered in
        # CONTACT_INFECTED_HOST_VECTOR).

    rates[RECOVER_HOST,:] = np.multiply( infected_hosts, recovery_rate_host )

    rates[RECOVER_VECTOR,:] = np.multiply(
        infected_vectors, recovery_rate_vector
        )

    rates[MUTATE_HOST,:] = np.multiply( infected_hosts, mutate_in_host )

    rates[MUTATE_VECTOR,:] = np.multiply( infected_vectors, mutate_in_vector )

    rates[RECOMBINE_HOST,:] = np.multiply( infected_hosts, recombine_in_host )

    rates[RECOMBINE_VECTOR,:] = np.multiply(
        infected_vectors, recombine_in_vector
        )

    rates[KILL_HOST,:] = np.multiply( infected_hosts, death_rate_host )

    rates[KILL_VECTOR,:] = np.multiply( infected_vectors, death_rate_vector )

    return rates
