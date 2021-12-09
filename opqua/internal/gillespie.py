
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

    MIGRATE_HOST = 0
    MIGRATE_VECTOR = 1
    POPULATION_CONTACT_HOST_HOST = 2
    POPULATION_CONTACT_HOST_VECTOR = 3
    POPULATION_CONTACT_VECTOR_HOST = 4
    CONTACT_HOST_HOST = 5
    CONTACT_HOST_VECTOR = 6
    CONTACT_VECTOR_HOST = 7
    RECOVER_HOST = 8
    RECOVER_VECTOR = 9
    MUTATE_HOST = 10
    MUTATE_VECTOR = 11
    RECOMBINE_HOST = 12
    RECOMBINE_VECTOR = 13
    KILL_HOST = 14
    KILL_VECTOR = 15
    DIE_HOST = 16
    DIE_VECTOR = 17
    BIRTH_HOST = 18
    BIRTH_VECTOR = 19

    EVENT_IDS = { # must match above
        0:'MIGRATE_HOST',
        1:'MIGRATE_VECTOR',
        2:'POPULATION_CONTACT_HOST_HOST',
        3:'POPULATION_CONTACT_HOST_VECTOR',
        4:'POPULATION_CONTACT_VECTOR_HOST',
        5:'CONTACT_HOST_HOST',
        6:'CONTACT_HOST_VECTOR',
        7:'CONTACT_VECTOR_HOST',
        8:'RECOVER_HOST',
        9:'RECOVER_VECTOR',
        10:'MUTATE_HOST',
        11:'MUTATE_VECTOR',
        12:'RECOMBINE_HOST',
        13:'RECOMBINE_VECTOR',
        14:'KILL_HOST',
        15:'KILL_VECTOR',
        16:'DIE_HOST',
        17:'DIE_VECTOR',
        18:'BIRTH_HOST',
        19:'BIRTH_VECTOR'
        }

    def __init__(self, model):
        """Create a new Gillespie simulation object.

        Arguments:
        model -- the model this simulation belongs to (Model)
        """

        super(Gillespie, self).__init__() # initialize as parent class object

        # Event IDs
        self.evt_IDs = [
            self.MIGRATE_HOST, self.MIGRATE_VECTOR,
            self.POPULATION_CONTACT_HOST_HOST,
            self.POPULATION_CONTACT_HOST_VECTOR,
            self.POPULATION_CONTACT_VECTOR_HOST,
            self.CONTACT_HOST_HOST, self.CONTACT_HOST_VECTOR,
            self.CONTACT_VECTOR_HOST,
            self.RECOVER_HOST, self.RECOVER_VECTOR,
            self.MUTATE_HOST, self.MUTATE_VECTOR,
            self.RECOMBINE_HOST, self.RECOMBINE_VECTOR,
            self.KILL_HOST, self.KILL_VECTOR,
            self.DIE_HOST, self.DIE_VECTOR,
            self.BIRTH_HOST, self.BIRTH_VECTOR
            ]
            # event IDs in specific order to be used

        self.total_population_contact_rate_host = 0
        self.total_population_contact_rate_vector = 0

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

        rates = np.zeros( [ len(self.evt_IDs), len(population_ids) ] )
            # rate array size of event space

        # Contact rates assume scaling area: large populations are equally
        # dense as small ones, so contact is constant with both host and
        # vector populations. If you don't want this to happen, modify each
        # population's base contact rate accordingly.

        # First calculate the total one-sided population contact rates for model
        self.total_population_contact_rate_host = np.sum([
            np.sum([ list(
                self.model.populations[id].neighbors_contact_hosts.values()
                ) ])
                * self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].RECEIVE_POPULATION_CONTACT
                    ].sum()
                / max( len( self.model.populations[id].hosts ), 1)
                for id in self.model.populations ])

        self.total_population_contact_rate_vector = np.sum([
            np.sum([ list(
                self.model.populations[id].neighbors_contact_vectors.values()
                ) ])
                * self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].RECEIVE_POPULATION_CONTACT
                    ].sum()
                / max( len( self.model.populations[id].vectors ), 1)
                for id in self.model.populations ])

        # Now for each population...
        for i,id in enumerate(population_ids):
            # Calculate the population's one-sided population contact rate
            self.model.populations[id].population_contact_rate_host = (
                np.sum([ list(
                    self.model.populations[id].neighbors_contact_hosts.values()
                    ) ])
                * self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].RECEIVE_POPULATION_CONTACT
                ].sum()
                / max( len( self.model.populations[id].hosts ), 1)
                )

            self.model.populations[id].population_contact_rate_vector = (
                np.sum([ list(
                    self.model.populations[id].neighbors_contact_vectors.values()
                    ) ])
                * self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].RECEIVE_POPULATION_CONTACT
                ].sum()
                / max( len( self.model.populations[id].vectors ), 1)
                )

            # Now the actual rates:
            rates[self.MIGRATE_HOST,i] = (
                self.model.populations[id].total_migration_rate_hosts
                * self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].MIGRATION
                    ].sum()
                )

            rates[self.MIGRATE_VECTOR,i] = (
                self.model.populations[id].total_migration_rate_vectors
                * self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].MIGRATION
                    ].sum()
                )

            rates[self.POPULATION_CONTACT_HOST_HOST,i] = (
                np.sum([ list(
                    self.model.populations[id].neighbors_contact_hosts.values()
                    ) ])
                * np.multiply(
                    self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].POPULATION_CONTACT
                        ],
                    self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].INFECTED
                        ]
                    ).sum()
                * ( self.total_population_contact_rate_host
                    -  self.model.populations[id].population_contact_rate_host )
                )

            rates[self.POPULATION_CONTACT_HOST_VECTOR,i] = (
                np.sum([ list(
                    self.model.populations[id].neighbors_contact_hosts.values()
                    ) ])
                * np.multiply(
                    self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].POPULATION_CONTACT
                        ],
                    self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].INFECTED
                        ]
                    ).sum()
                * ( self.total_population_contact_rate_vector
                    -self.model.populations[id].population_contact_rate_vector )
                )

            rates[self.POPULATION_CONTACT_VECTOR_HOST,i] = (
                np.sum([ list(
                    self.model.populations[id].neighbors_contact_vectors.values()
                    ) ])
                * np.multiply(
                    self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].POPULATION_CONTACT
                        ],
                    self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].INFECTED
                        ]
                    ).sum()
                * ( self.total_population_contact_rate_host
                    - self.model.populations[id].population_contact_rate_host )
                )

            rates[self.CONTACT_HOST_HOST,i] = (
                self.model.populations[id].contact_rate_host_host
                * np.multiply(
                    self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].CONTACT
                        ],
                    self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].INFECTED
                        ]
                    ).sum()
                * self.model.populations[id].transmission_efficiency_host_host
                * np.sum( self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].RECEIVE_CONTACT
                    ] )
                / max( len( self.model.populations[id].hosts ), 1)
                )

            rates[self.CONTACT_HOST_VECTOR,i] = (
                self.model.populations[id].contact_rate_host_vector
                * np.multiply(
                    self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].CONTACT
                        ],
                    self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].INFECTED
                        ]
                    ).sum()
                * self.model.populations[id].transmission_efficiency_host_vector
                * np.sum( self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].RECEIVE_CONTACT
                    ] )
                / max( len( self.model.populations[id].hosts ), 1)
                )

            rates[self.CONTACT_VECTOR_HOST,i] = (
                self.model.populations[id].contact_rate_host_vector
                * np.multiply(
                    self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].CONTACT
                        ],
                    self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].INFECTED
                        ]
                    ).sum()
                * self.model.populations[id].transmission_efficiency_vector_host
                * np.sum( self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].RECEIVE_CONTACT
                    ])
                / max( len( self.model.populations[id].hosts ), 1)
                )

            rates[self.RECOVER_HOST,i] = (
                self.model.populations[id].recovery_rate_host
                * self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].RECOVERY
                    ].sum()
                )

            rates[self.RECOVER_VECTOR,i] = (
                self.model.populations[id].recovery_rate_vector
                * self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].RECOVERY
                    ].sum()
                )

            rates[self.MUTATE_HOST,i] = (
                self.model.populations[id].mutate_in_host
                * self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].MUTATION
                    ].sum()
                )

            rates[self.MUTATE_VECTOR,i] = (
                self.model.populations[id].mutate_in_vector
                * self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].MUTATION
                    ].sum()
                )

            rates[self.RECOMBINE_HOST,i] = (
                self.model.populations[id].recombine_in_host
                * self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].RECOMBINATION
                    ].sum()
                )

            rates[self.RECOMBINE_VECTOR,i] = (
                self.model.populations[id].recombine_in_vector
                * self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].RECOMBINATION
                    ].sum()
                )

            rates[self.KILL_HOST,i] = (
                self.model.populations[id].lethality_rate_host
                * self.model.populations[id].recovery_rate_host
                * self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].LETHALITY
                    ].sum()
                )

            rates[self.KILL_VECTOR,i] = (
                self.model.populations[id].lethality_rate_vector
                * self.model.populations[id].recovery_rate_vector
                * self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].LETHALITY
                    ].sum()
                )

            rates[self.DIE_HOST,i] = (
                self.model.populations[id].death_rate_host
                * len(self.model.populations[id].hosts)
                )

            rates[self.DIE_VECTOR,i] = (
                self.model.populations[id].death_rate_vector
                * len(self.model.populations[id].vectors)
                )

            rates[self.BIRTH_HOST,i] = (
                self.model.populations[id].birth_rate_host
                * self.model.populations[id].coefficients_hosts[
                    :, self.model.populations[id].NATALITY
                    ].sum()
                )

            rates[self.BIRTH_VECTOR,i] = (
                self.model.populations[id].birth_rate_vector
                * self.model.populations[id].coefficients_vectors[
                    :, self.model.populations[id].NATALITY
                    ].sum()
                )

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

        if act == self.MIGRATE_HOST:
            rand = rand * pop.total_migration_rate_hosts
            r_cum = 0
            for neighbor in pop.neighbors_hosts:
                r_cum += pop.neighbors_hosts[neighbor]
                if r_cum > rand:
                    pop.migrate(neighbor,1,0, rand=(
                        ( rand - r_cum + pop.neighbors_hosts[neighbor] )
                        / pop.neighbors_hosts[neighbor] ) )
                    changed = True
                    break

        elif act == self.MIGRATE_VECTOR:
            rand = rand * pop.total_migration_rate_vectors
            r_cum = 0
            for neighbor in pop.neighbors_vectors:
                r_cum += pop.neighbors_vectors[neighbor]
                if r_cum > rand:
                    pop.migrate(neighbor,0,1, rand=(
                        ( rand - r_cum + pop.neighbors_vectors[neighbor] )
                        / pop.neighbors_vectors[neighbor] ) )
                    changed = True
                    break

        elif act == self.POPULATION_CONTACT_HOST_HOST:
            rand = rand * np.sum([
                pop.neighbors_contact_hosts[neighbor]
                * neighbor.neighbors_contact_hosts[pop]
                * neighbor.coefficients_hosts[
                    :, neighbor.RECEIVE_POPULATION_CONTACT
                    ].sum()
                / max( len( neighbor.hosts ), 1)
                for neighbor in pop.neighbors_contact_hosts ])
            r_cum = 0
            for neighbor in pop.neighbors_contact_hosts:
                neighbor_r = (
                    pop.neighbors_contact_hosts[neighbor]
                    * neighbor.neighbors_contact_hosts[pop]
                    * neighbor.coefficients_hosts[
                        :, neighbor.RECEIVE_POPULATION_CONTACT
                        ].sum()
                    / max( len( neighbor.hosts ), 1)
                    )
                r_cum += neighbor_r
                if r_cum > rand:
                    changed = pop.populationContact(
                        neighbor, ( rand - r_cum + neighbor_r ) / neighbor_r,
                        host_origin=True, host_target=True
                        )
                    break

        elif act == self.POPULATION_CONTACT_HOST_VECTOR:
            rand = rand * np.sum([
                pop.neighbors_contact_hosts[neighbor]
                * neighbor.neighbors_contact_vectors[pop]
                * neighbor.coefficients_vectors[
                    :, neighbor.RECEIVE_POPULATION_CONTACT
                    ].sum()
                / max( len( neighbor.vectors ), 1)
                for neighbor in pop.neighbors_contact_hosts ])
            r_cum = 0
            for neighbor in pop.neighbors_contact_hosts:
                neighbor_r = (
                    pop.neighbors_contact_hosts[neighbor]
                    * neighbor.neighbors_contact_vectors[pop]
                    * neighbor.coefficients_vectors[
                        :, neighbor.RECEIVE_POPULATION_CONTACT
                        ].sum()
                    / max( len( neighbor.vectors ), 1)
                    )
                r_cum += neighbor_r
                if r_cum > rand:
                    changed = pop.populationContact(
                        neighbor, ( rand - r_cum + neighbor_r ) / neighbor_r,
                        host_origin=True, host_target=False
                        )
                    break

        elif act == self.POPULATION_CONTACT_VECTOR_HOST:
            rand = rand * np.sum([
                pop.neighbors_contact_vectors[neighbor]
                * neighbor.neighbors_contact_hosts[pop]
                * neighbor.coefficients_hosts[
                    :, neighbor.RECEIVE_POPULATION_CONTACT
                    ].sum()
                / max( len( neighbor.hosts ), 1)
                for neighbor in pop.neighbors_contact_vectors ])
            r_cum = 0
            for neighbor in pop.neighbors_contact_vectors:
                neighbor_r = (
                    pop.neighbors_contact_vectors[neighbor]
                    * neighbor.neighbors_contact_hosts[pop]
                    * neighbor.coefficients_hosts[
                        :, neighbor.POPULATION_CONTACT
                        ].sum()
                    / max( len( neighbor.hosts ), 1)
                    )
                r_cum += neighbor_r
                if r_cum > rand:
                    changed = pop.populationContact(
                        neighbor, ( rand - r_cum + neighbor_r ) / neighbor_r,
                        host_origin=False, host_target=True
                        )
                    break

        elif act == self.CONTACT_HOST_HOST:
            changed = pop.contactHostHost(rand)

        elif act == self.CONTACT_HOST_VECTOR:
            changed = pop.contactHostVector(rand)

        elif act == self.CONTACT_VECTOR_HOST:
            changed = pop.contactVectorHost(rand)

        elif act == self.RECOVER_HOST:
            pop.recoverHost(rand)
            changed = True

        elif act == self.RECOVER_VECTOR:
            pop.recoverVector(rand)
            changed = True

        elif act == self.MUTATE_HOST:
            pop.mutateHost(rand)
            changed = True

        elif act == self.MUTATE_VECTOR:
            pop.mutateVector(rand)
            changed = True

        elif act == self.RECOMBINE_HOST:
            pop.recombineHost(rand)
            changed = True

        elif act == self.RECOMBINE_VECTOR:
            pop.recombineVector(rand)
            changed = True

        elif act == self.KILL_HOST:
            pop.killHost(rand)
            changed = True

        elif act == self.KILL_VECTOR:
            pop.killVector(rand)
            changed = True

        elif act == self.DIE_HOST:
            pop.dieHost(rand)
            changed = True

        elif act == self.DIE_VECTOR:
            pop.dieVector(rand)
            changed = True

        elif act == self.BIRTH_HOST:
            pop.birthHost(rand)
            changed = True

        elif act == self.BIRTH_VECTOR:
            pop.birthVector(rand)
            changed = True

        self.model.global_trackers['num_events'][self.EVENT_IDS[act]] += 1

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
        self.model.t_var = t0 # keeps track of time
        history = { 0: self.model.copyState(
            host_sampling=host_sampling,
            vector_sampling=vector_sampling
            ) }
        intervention_tracker = 0
            # keeps track of what the next intervention should be
        self.model.interventions = sorted(
            self.model.interventions, key=lambda i: i.time
            )

        print_counter = 0 # only used to track when to print
        sampling_counter = 0 # used to track when to save a snapshot

        while self.model.t_var < tf: # repeat until t reaches end of timecourse
            population_ids = list( self.model.populations.keys() )
            r = self.getRates(population_ids) # get event rates in this state
            r_tot = np.sum(r) # sum of all rates

            # Time handling
            if r_tot > 0:
                dt = np.random.exponential( 1/r_tot ) # time until next event

                if (intervention_tracker < len(self.model.interventions)
                    and self.model.t_var + dt
                    >= self.model.interventions[intervention_tracker].time):
                    # if there are any interventions left and if it is time
                    # to make one,
                    while ( intervention_tracker < len(self.model.interventions)
                        and (self.model.t_var + dt
                        >= self.model.interventions[intervention_tracker].time
                        or r_tot == 0) and self.model.t_var < tf ):
                            # carry out all interventions at this time point,
                            # and additional timepoints if no events will happen
                        self.model.t_var = self.model.interventions[
                            intervention_tracker
                            ].time
                        self.model.interventions[
                            intervention_tracker
                            ].doIntervention()
                        intervention_tracker += 1 # advance the tracker

                        # save snapshot at this timepoint
                        sampling_counter = 0
                        history[self.model.t_var] = self.model.copyState()

                        # now recalculate rates
                        population_ids = list( self.model.populations.keys() )
                        r = self.getRates(population_ids)
                            # get event rates in this state
                        r_tot = np.sum(r) # sum of all rates

                    if r_tot > 0: # if no more events happening,
                        dt = np.random.exponential( 1/r_tot )
                            # time until next event
                    else:
                        self.model.t_var = tf # go to end

                self.model.t_var += dt # add time step to main timer

                # Event handling
                if self.model.t_var < tf: # if still within max time
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
                                        + str(self.model.t_var) + ', event: '
                                        + self.EVENT_IDS[e]
                                        )

                                changed = self.doAction(
                                    e, self.model.populations[
                                        population_ids[p]
                                        ], ( u - r_cum + r[e,p] ) / r[e,p]
                                    ) # do corresponding action,
                                      # feed in renormalized random number

                                if changed: # if model state changed
                                    # update custom condition trackers
                                    for condition in self.model.custom_condition_trackers:
                                        if self.model.custom_condition_trackers[condition](self.model):
                                            self.model.global_trackers['custom_conditions'][condition].append(self.model.t_var)

                                    if time_sampling >= 0:
                                            # if state changed and saving history,
                                            # saves history at correct intervals
                                        sampling_counter += 1
                                        if sampling_counter > time_sampling:
                                            sampling_counter = 0
                                            history[self.model.t_var] = self.model.copyState(
                                                host_sampling=host_sampling,
                                                vector_sampling=vector_sampling
                                                )

                                break # exit event loop

                        else: # if the inner loop wasn't broken,
                            continue # continue outer loop

                        break # otherwise, break outer loop
            else: # if no events happening,
                if intervention_tracker < len(self.model.interventions):
                        # if still not done with interventions,
                    while (intervention_tracker < len(self.model.interventions)
                        and self.model.t_var
                        <= self.model.interventions[intervention_tracker].time
                        and self.model.t_var < tf):
                            # carry out all interventions at this time point
                        self.model.t_var = self.model.interventions[
                            intervention_tracker
                            ].time
                        self.model.interventions[
                            intervention_tracker
                            ].doIntervention()
                        intervention_tracker += 1 # advance the tracker
                else:
                    self.model.t_var = tf

        print( 'Simulating time: ' + str(self.model.t_var), 'END')
        history[tf] = self.model.copyState(
            host_sampling=host_sampling,
            vector_sampling=vector_sampling
            )
        history[tf].history = None
        history[tf].global_trackers = cp.deepcopy( self.model.global_trackers )

        return history
