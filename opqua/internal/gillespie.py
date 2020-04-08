
"""Contains class Population."""

import copy as cp
import pandas as pd
import numpy as np
import joblib as jl

from opqua.model import *

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

    def __init__(self, model):
        """Create a new Gillespie simulation object.

        Arguments:
        model -- the model this simulation belongs to (Model)
        """

        super(Gillespie, self).__init__() # initialize as parent class object

        # Event IDs
        self.evt_IDs = [
            self.MIGRATE, self.CONTACT_INFECTED_HOST_ANY_HOST,
            self.CONTACT_INFECTED_HOST_ANY_VECTOR,
            self.CONTACT_HEALTHY_HOST_INFECTED_VECTOR,
            self.RECOVER_HOST, self.RECOVER_VECTOR,
            self.MUTATE_HOST, self.MUTATE_VECTOR,
            self.RECOMBINE_HOST, self.RECOMBINE_VECTOR,
            self.KILL_HOST, self.KILL_VECTOR
            ]
            # event IDs in specific order to be used

        self.model = model

    def getRates(self,population_ids):
        """Calculates event rates according to current system state.

        Arguments:
        population_ids -- list with ids for every population in the model
            (list of Strings)

        Returns:
        dictionary with event ID constants as keys and rates as values. Includes
            total rate under 'tot' key.
        """

        rates = np.zeros( [ len(self.evt_IDs), len(population_ids) ] )
            # rate array size of event space

        rates[self.MIGRATE,:] = np.array( [
            len( self.model.populations[id].hosts )
            * self.model.populations[id].total_migration_rate
            for id in population_ids
            ] )

        rates[self.CONTACT_INFECTED_HOST_ANY_HOST,:] = np.array( [
            ( len( self.model.populations[id].infected_hosts )
            / np.maximum( len( self.model.populations[id].hosts ), 1 ) )
            * self.model.populations[id].contact_rate_host_host
            for id,p in self.model.populations.items()
            ] )
            # contact rates assumes scaling area--large populations are equally
            # dense as small ones, so contact is constant with both host and
            # vector populations. If you don't want this to happen, modify the
            # population's contact rate accordingly.
            # Examines contacts between infected hosts and all hosts

        rates[self.CONTACT_INFECTED_HOST_ANY_VECTOR,:] = np.array( [
            ( len( self.model.populations[id].infected_hosts )
            / np.maximum( len( self.model.populations[id].hosts ), 1 ) )
            * self.model.populations[id].contact_rate_host_vector
            for id,p in self.model.populations.items()
            ] )
            # contact rates assumes scaling area--large populations are equally
            # dense as small ones, so contact is constant with both host and
            # vector populations. If you don't want this to happen, modify the
            # population's contact rate accordingly.
            # Examines contacts between infected hosts and all vectors

        rates[self.CONTACT_HEALTHY_HOST_INFECTED_VECTOR,:] = np.array( [
            ( len( self.model.populations[id].healthy_hosts )
            / np.maximum( len( self.model.populations[id].hosts ), 1 ) )
            * ( len( self.model.populations[id].infected_vectors )
            / np.maximum( len( self.model.populations[id].vectors ), 1 ) )
            * self.model.populations[id].contact_rate_host_vector
            for id,p in self.model.populations.items()
            ] )
            # contact rates assumes scaling area--large populations are equally
            # dense as small ones, so contact is constant with both host and
            # vector populations. If you don't want this to happen, modify the
            # population's contact rate accordingly.
            # Examines contacts between healthy hosts and infected vectors
            # (infected-infected contacts are considered in
            # CONTACT_INFECTED_HOST_VECTOR).

        rates[self.RECOVER_HOST,:] = np.array( [
            len( self.model.populations[id].infected_hosts )
            * self.model.populations[id].recovery_rate_host
            for id,p in self.model.populations.items()
            ] )

        rates[self.RECOVER_VECTOR,:] = np.array( [
            len( self.model.populations[id].infected_vectors )
            * self.model.populations[id].recovery_rate_vector
            for id,p in self.model.populations.items()
            ] )

        rates[self.MUTATE_HOST,:] = np.array( [
            len( self.model.populations[id].infected_hosts )
            * self.model.populations[id].mutate_in_host
            for id,p in self.model.populations.items()
            ] )

        rates[self.MUTATE_VECTOR,:] = np.array( [
            len( self.model.populations[id].infected_vectors )
            * self.model.populations[id].mutate_in_vector
            for id,p in self.model.populations.items()
            ] )

        rates[self.RECOMBINE_HOST,:] = np.array( [
            len( self.model.populations[id].infected_hosts )
            * self.model.populations[id].recombine_in_host
            for id,p in self.model.populations.items()
            ] )

        rates[self.RECOMBINE_VECTOR,:] = np.array( [
            len( self.model.populations[id].infected_vectors )
            * self.model.populations[id].recombine_in_vector
            for id,p in self.model.populations.items()
            ] )

        rates[self.KILL_HOST,:] = np.array( [
            len( self.model.populations[id].infected_hosts )
            * self.model.populations[id].death_rate_host
            for id,p in self.model.populations.items()
            ] )

        rates[self.KILL_VECTOR,:] = np.array( [
            len( self.model.populations[id].infected_vectors )
            * self.model.populations[id].death_rate_vector
            for id,p in self.model.populations.items()
            ] )

        return rates

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

        if act == self.MIGRATE:
            rand = rand * pop.total_migration_rate
            r_cum = 0
            for neighbor in pop.neighbors:
                r_cum += pop.neighbors[neighbor]
                if r_cum > rand:
                    pop.migrate(neighbor,1,0)
                    changed = True
                    break



        elif act == self.CONTACT_INFECTED_HOST_ANY_HOST:
            rand = rand * len(pop.infected_hosts)
            infected_host = int( np.floor(rand) )
            other_host = int(
                np.floor( ( rand - infected_host ) * len(pop.hosts) )
                )
            changed = pop.contactInfectedHostAnyHost(infected_host,other_host)

        elif act == self.CONTACT_INFECTED_HOST_ANY_VECTOR:
            rand = rand * len(pop.infected_hosts)
            infected_host = int( np.floor(rand) )
            vector = int(
                np.floor( ( rand - infected_host ) * len(pop.vectors) )
                )
            changed = pop.contactInfectedHostAnyVector(infected_host,vector)

        elif act == self.CONTACT_HEALTHY_HOST_INFECTED_VECTOR:
            rand = rand * len(pop.healthy_hosts)
            healthy_host = int( np.floor(rand) )
            infected_vector = int(
                np.floor( ( rand - healthy_host ) * len(pop.infected_vectors) )
                )
            changed = pop.contactHealthyHostInfectedVector(
                healthy_host,infected_vector
                )

        elif act == self.RECOVER_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            pop.recoverHost(host)
            changed = True

        elif act == self.RECOVER_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.recoverVector(vector)
            changed = True

        elif act == self.MUTATE_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            pop.mutateHost(host)
            changed = True

        elif act == self.MUTATE_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.mutateVector(vector)
            changed = True

        elif act == self.RECOMBINE_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            pop.recombineHost(host)
            changed = True

        elif act == self.RECOMBINE_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.recombineVector(vector)
            changed = True

        elif act == self.KILL_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            pop.killHost(host)
            changed = True

        elif act == self.KILL_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.killVector(vector)
            changed = True

        return changed

    def run(self,t0,tf):

        """Simulate model for a specified time between two time points.

        Simulates a time series using the Gillespie algorithm.

        Arguments:
        t0 -- initial time point to start simulation at (number)
        tf -- initial time point to end simulation at (number)

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
                                if changed:
                                    history[t_var] = cp.deepcopy(self.model)

                                break # exit event loop

                        else: # if the inner loop wasn't broken,
                            continue # continue outer loop

                        break # otherwise, break outer loop
            else: # if at end of time course,
                if intervention_tracker < len(self.model.interventions):
                    print( 'Simulating time: ' + str(t_var), e)
                    self.model.interventions[
                        intervention_tracker
                        ].doIntervention()
                    t_var = self.model.interventions[intervention_tracker].time
                    intervention_tracker += 1 # advance the tracker
                else:
                    print( 'Simulating time: ' + str(t_var), e)
                    t_var = tf

        history[tf] = cp.deepcopy(self.model)
        history[tf].history = None

        return history
