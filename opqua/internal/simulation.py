
"""Contains class Population."""

import copy as cp
import pandas as pd
import numpy as np
import joblib as jl

from opqua.model import *

class Simulation(object):
    """Class contains methods for simulating model.

    Class defines a model's events and methods for changing system state
    according to the possible events and simulating a timecourse using the
    exact Gillespie algorithm or an approximated, mixed Euler-Gillespie model
    (a variation of tau-leaping).

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
        MIGRATE_HOST:'MIGRATE_HOST',
        MIGRATE_VECTOR:'MIGRATE_VECTOR',
        POPULATION_CONTACT_HOST_HOST:'POPULATION_CONTACT_HOST_HOST',
        POPULATION_CONTACT_HOST_VECTOR:'POPULATION_CONTACT_HOST_VECTOR',
        POPULATION_CONTACT_VECTOR_HOST:'POPULATION_CONTACT_VECTOR_HOST',
        CONTACT_HOST_HOST:'CONTACT_HOST_HOST',
        CONTACT_HOST_VECTOR:'CONTACT_HOST_VECTOR',
        CONTACT_VECTOR_HOST:'CONTACT_VECTOR_HOST',
        RECOVER_HOST:'RECOVER_HOST',
        RECOVER_VECTOR:'RECOVER_VECTOR',
        MUTATE_HOST:'MUTATE_HOST',
        MUTATE_VECTOR:'MUTATE_VECTOR',
        RECOMBINE_HOST:'RECOMBINE_HOST',
        RECOMBINE_VECTOR:'RECOMBINE_VECTOR',
        KILL_HOST:'KILL_HOST',
        KILL_VECTOR:'KILL_VECTOR',
        DIE_HOST:'DIE_HOST',
        DIE_VECTOR:'DIE_VECTOR',
        BIRTH_HOST:'BIRTH_HOST',
        BIRTH_VECTOR:'BIRTH_VECTOR'
        }

    def __init__(self, model):
        """Create a new simulation object.

        Arguments:
        model -- the model this simulation belongs to (Model)
        """

        super(Simulation, self).__init__() # initialize as parent class object

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

        self.evt_possible = np.repeat( 1, len(self.evt_IDs) )
            # True/False whether or not each event might possibly happen in this
            # simulation: we will fill it out in the next lines

        interventions = [ i.intervention for i in model.interventions ]

        if ( np.sum([
                    len(pop.vectors) for i,pop in model.populations.items()
                    ]) == 0 and 'addVectors' not in interventions ):
            self.evt_possible[self.MIGRATE_VECTOR] = 0
            self.evt_possible[self.POPULATION_CONTACT_HOST_VECTOR] = 0
            self.evt_possible[self.POPULATION_CONTACT_VECTOR_HOST] = 0
            self.evt_possible[self.CONTACT_HOST_VECTOR] = 0
            self.evt_possible[self.RECOVER_VECTOR] = 0
            self.evt_possible[self.MUTATE_VECTOR] = 0
            self.evt_possible[self.RECOMBINE_VECTOR] = 0
            self.evt_possible[self.KILL_VECTOR] = 0
            self.evt_possible[self.DIE_VECTOR] = 0
            self.evt_possible[self.BIRTH_VECTOR] = 0
        if ( len(model.populations) < 2 ):
            self.evt_possible[self.MIGRATE_HOST] = 0
            self.evt_possible[self.MIGRATE_VECTOR] = 0
            self.evt_possible[self.POPULATION_CONTACT_HOST_HOST] = 0
            self.evt_possible[self.POPULATION_CONTACT_HOST_VECTOR] = 0
            self.evt_possible[self.POPULATION_CONTACT_VECTOR_HOST] = 0
        if ( np.sum([
                s.contact_rate_host_host for i,s in model.setups.items()
                ]) == 0 ):
            self.evt_possible[self.CONTACT_HOST_HOST] = 0
        if ( np.sum([
                s.contact_rate_host_vector for i,s in model.setups.items()
                ]) == 0 ):
            self.evt_possible[self.CONTACT_HOST_VECTOR] = 0
            self.evt_possible[self.CONTACT_VECTOR_HOST] = 0
        if ( np.sum([
                s.recovery_rate_host for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.RECOVER_HOST] = 0
        if ( np.sum([
                s.recovery_rate_vector for i,s in model.setups.items()
                ]) == 0 ):
            self.evt_possible[self.RECOVER_VECTOR] = 0
        if ( np.sum([
                s.mutate_in_host for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.MUTATE_HOST] = 0
        if ( np.sum([
                s.mutate_in_vector for i,s in model.setups.items()  ]) == 0 ):
            self.evt_possible[self.MUTATE_VECTOR] = 0
        if ( np.sum([
                s.recombine_in_host for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.RECOMBINE_HOST] = 0
        if ( np.sum([
                s.recombine_in_vector for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.RECOMBINE_VECTOR] = 0
        if ( np.sum([
                s.mortality_rate_host for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.KILL_HOST] = 0
        if ( np.sum([
                s.mortality_rate_vector for i,s in model.setups.items()
                ]) == 0 ):
            self.evt_possible[self.KILL_VECTOR] = 0
        if ( np.sum([
                s.death_rate_host for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.DIE_HOST] = 0
        if ( np.sum([
                s.death_rate_vector for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.DIE_VECTOR] = 0
        if ( np.sum([
                s.birth_rate_host for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.BIRTH_HOST] = 0
        if ( np.sum([
                s.birth_rate_vector for i,s in model.setups.items() ]) == 0 ):
            self.evt_possible[self.BIRTH_VECTOR] = 0

        self.total_population_contact_rate_host = 0
        self.total_population_contact_rate_vector = 0

        self.model = model
        self.random = model.random # random number generator

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

        # Now for each population...
        for i,id in enumerate(population_ids):
            # Calculate the rates:
            if self.evt_possible[self.MIGRATE_HOST]:
                rates[self.MIGRATE_HOST,i] = (
                    self.model.populations[id].total_migration_rate_hosts
                    * self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].MIGRATION
                        ].sum()
                    )

            if self.evt_possible[self.MIGRATE_VECTOR]:
                rates[self.MIGRATE_VECTOR,i] = (
                    self.model.populations[id].total_migration_rate_vectors
                    * self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].MIGRATION
                        ].sum()
                    )

            if self.evt_possible[self.POPULATION_CONTACT_HOST_HOST]:
                rates[self.POPULATION_CONTACT_HOST_HOST,i] = (
                    np.sum([ list(
                        self.model.populations[id].neighbors_contact_hosts_hosts.values()
                        ) ])
                    * np.multiply(
                        self.model.populations[id].coefficients_hosts[
                            :, self.model.populations[id].POPULATION_CONTACT
                            ],
                        self.model.populations[id].coefficients_hosts[
                            :, self.model.populations[id].INFECTED
                            ]
                        ).sum() / max( len( self.model.populations[id].hosts ), 1)
                    * np.sum([
                        neighbor.neighbors_contact_hosts_hosts[self.model.populations[id]]
                        * neighbor.coefficients_hosts[
                            :, neighbor.RECEIVE_POPULATION_CONTACT
                            ].sum()
                        for neighbor in self.model.populations[id].neighbors_contact_hosts_hosts
                        ])
                    )

            if self.evt_possible[self.POPULATION_CONTACT_HOST_VECTOR]:
                rates[self.POPULATION_CONTACT_HOST_VECTOR,i] = (
                    np.sum([ list(
                        self.model.populations[id].neighbors_contact_hosts_vectors.values()
                        ) ])
                    * np.multiply(
                        self.model.populations[id].coefficients_hosts[
                            :, self.model.populations[id].POPULATION_CONTACT
                            ],
                        self.model.populations[id].coefficients_hosts[
                            :, self.model.populations[id].INFECTED
                            ]
                        ).sum() / max( len( self.model.populations[id].hosts ), 1)
                    * np.sum([
                        neighbor.neighbors_contact_vectors_hosts[self.model.populations[id]]
                        * neighbor.coefficients_vectors[
                            :, neighbor.RECEIVE_POPULATION_CONTACT
                            ].sum()
                        for neighbor in self.model.populations[id].neighbors_contact_hosts_vectors
                        ])
                    )

            if self.evt_possible[self.POPULATION_CONTACT_VECTOR_HOST]:
                rates[self.POPULATION_CONTACT_VECTOR_HOST,i] = (
                    np.sum([ list(
                        self.model.populations[id].neighbors_contact_vectors_hosts.values()
                        ) ])
                    * np.multiply(
                        self.model.populations[id].coefficients_vectors[
                            :, self.model.populations[id].POPULATION_CONTACT
                            ],
                        self.model.populations[id].coefficients_vectors[
                            :, self.model.populations[id].INFECTED
                            ]
                        ).sum()
                    * np.sum([
                        neighbor.neighbors_contact_hosts_vectors[self.model.populations[id]]
                        * neighbor.coefficients_hosts[
                            :, neighbor.RECEIVE_POPULATION_CONTACT
                            ].sum() / max( len( neighbor.hosts ), 1)
                        for neighbor in self.model.populations[id].neighbors_contact_vectors_hosts
                        ])
                    )

            if self.evt_possible[self.CONTACT_HOST_HOST]:
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

            if self.evt_possible[self.CONTACT_HOST_VECTOR]:
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

            if self.evt_possible[self.CONTACT_VECTOR_HOST]:
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

            if self.evt_possible[self.RECOVER_HOST]:
                rates[self.RECOVER_HOST,i] = (
                    self.model.populations[id].recovery_rate_host
                    * self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].RECOVERY
                        ].sum()
                    )

            if self.evt_possible[self.RECOVER_VECTOR]:
                rates[self.RECOVER_VECTOR,i] = (
                    self.model.populations[id].recovery_rate_vector
                    * self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].RECOVERY
                        ].sum()
                    )

            if self.evt_possible[self.MUTATE_HOST]:
                rates[self.MUTATE_HOST,i] = (
                    self.model.populations[id].mutate_in_host
                    * self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].MUTATION
                        ].sum()
                    )

            if self.evt_possible[self.MUTATE_VECTOR]:
                rates[self.MUTATE_VECTOR,i] = (
                    self.model.populations[id].mutate_in_vector
                    * self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].MUTATION
                        ].sum()
                    )

            if self.evt_possible[self.RECOMBINE_HOST]:
                rates[self.RECOMBINE_HOST,i] = (
                    self.model.populations[id].recombine_in_host
                    * self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].RECOMBINATION
                        ].sum()
                    )

            if self.evt_possible[self.RECOMBINE_VECTOR]:
                rates[self.RECOMBINE_VECTOR,i] = (
                    self.model.populations[id].recombine_in_vector
                    * self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].RECOMBINATION
                        ].sum()
                    )

            if self.evt_possible[self.KILL_HOST]:
                rates[self.KILL_HOST,i] = (
                    self.model.populations[id].mortality_rate_host
                    * self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].LETHALITY
                        ].sum()
                    )

            if self.evt_possible[self.KILL_VECTOR]:
                rates[self.KILL_VECTOR,i] = (
                    self.model.populations[id].mortality_rate_vector
                    * self.model.populations[id].coefficients_vectors[
                        :, self.model.populations[id].LETHALITY
                        ].sum()
                    )

            if self.evt_possible[self.DIE_HOST]:
                rates[self.DIE_HOST,i] = (
                    self.model.populations[id].death_rate_host
                    * len(self.model.populations[id].hosts)
                    )

            if self.evt_possible[self.DIE_VECTOR]:
                rates[self.DIE_VECTOR,i] = (
                    self.model.populations[id].death_rate_vector
                    * len(self.model.populations[id].vectors)
                    )

            if self.evt_possible[self.BIRTH_HOST]:
                rates[self.BIRTH_HOST,i] = (
                    self.model.populations[id].birth_rate_host
                    * self.model.populations[id].coefficients_hosts[
                        :, self.model.populations[id].NATALITY
                        ].sum()
                    )

            if self.evt_possible[self.BIRTH_VECTOR]:
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
            neighbors = list( pop.neighbors_contact_hosts_hosts.keys() )
            neighbor_rates = [
                pop.neighbors_contact_hosts_hosts[neighbor]
                * neighbor.neighbors_contact_hosts_hosts[pop]
                * neighbor.coefficients_hosts[
                    :, neighbor.RECEIVE_POPULATION_CONTACT
                    ].sum()
                / max( len( pop.hosts ), 1) # technically unnecessary
                for neighbor in neighbors
                ]
            rand = rand * np.sum(neighbor_rates)
            r_cum = 0
            for neighbor_i,neighbor_r in enumerate(neighbor_rates):
                r_cum += neighbor_r
                if r_cum > rand:
                    changed = pop.populationContact(
                        neighbors[neighbor_i],
                        ( rand - r_cum + neighbor_r ) / neighbor_r,
                        host_origin=True, host_target=True
                        )
                    break

        elif act == self.POPULATION_CONTACT_HOST_VECTOR:
            neighbors = list( pop.neighbors_contact_hosts_vectors.keys() )
            neighbor_rates = [
                pop.neighbors_contact_hosts_vectors[neighbor]
                * neighbor.neighbors_contact_vectors_hosts[pop]
                * neighbor.coefficients_vectors[
                    :, neighbor.RECEIVE_POPULATION_CONTACT
                    ].sum()
                / max( len( pop.hosts ), 1) # technically unnecessary
                for neighbor in neighbors
                ]
            rand = rand * np.sum(neighbor_rates)
            r_cum = 0
            for neighbor_i,neighbor_r in enumerate(neighbor_rates):
                r_cum += neighbor_r
                if r_cum > rand:
                    changed = pop.populationContact(
                        neighbors[neighbor_i],
                        ( rand - r_cum + neighbor_r ) / neighbor_r,
                        host_origin=True, host_target=False
                        )
                    break

        elif act == self.POPULATION_CONTACT_VECTOR_HOST:
            neighbors = list( pop.neighbors_contact_vectors_hosts.keys() )
            neighbor_rates = [
                pop.neighbors_contact_vectors_hosts[neighbor]
                * neighbor.neighbors_contact_hosts_vectors[pop]
                * neighbor.coefficients_hosts[
                    :, neighbor.RECEIVE_POPULATION_CONTACT
                    ].sum()
                / max( len( neighbor.hosts ), 1) # necessary!
                for neighbor in neighbors
                ]
            rand = rand * np.sum(neighbor_rates)
            r_cum = 0
            for neighbor_i,neighbor_r in enumerate(neighbor_rates):
                r_cum += neighbor_r
                if r_cum > rand:
                    changed = pop.populationContact(
                        neighbors[neighbor_i],
                        ( rand - r_cum + neighbor_r ) / neighbor_r,
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

    def run(self,t0,tf,method='exact',dt_leap=None,dt_thre=None,
            time_sampling=0,host_sampling=0,vector_sampling=0,
            skip_uninfected=False,print_every_n_events=1000):
        """Simulate model for a specified time between two time points.

        Wrapper for different simulation algorithms.

        Arguments:
        t0 -- initial time point to start simulation at (number)
        tf -- initial time point to end simulation at (number)
        method -- algorithm to be used; default is approximated solver
            ('approximated' or 'exact'; default 'exact')
        dt_leap -- time leap size used to simulate bursts; if None, set to
            minimum growth threshold time across all populations (number,
            default None)
        dt_thre -- time threshold below which bursts are used; if None, set to
            dt_leap (number, default None)
        time_sampling -- how many events to skip before saving a snapshot of the
            system state (saves all by default), if <0, saves only final state
            (int, default 0)
        host_sampling -- how many hosts to skip before saving one in a snapshot
            of the system state (saves all by default) (int, default 0)
        vector_sampling -- how many vectors to skip before saving one in a
            snapshot of the system state (saves all by default) (int, default 0)
        skip_uninfected -- whether to save only infected hosts/vectors and
            record the number of uninfected host/vectors instead (Boolean,
            default False)
        print_every_n_events -- number of events a message is printed to console
            (int>0, default 1000)

        Returns:
        dictionary containing model state history, with keys=times and
            values=Model objects with model snapshot at that time point
        """

        if method=='approximated':
            return self.runApproximated(t0,tf,dt_leap,dt_thre,time_sampling,
                host_sampling,vector_sampling,skip_uninfected,
                print_every_n_events)
        elif method=='exact':
            return self.runExactGillespie(t0,tf,time_sampling,host_sampling,
                vector_sampling,skip_uninfected,print_every_n_events)


    def runExactGillespie(self,t0,tf,time_sampling=0,host_sampling=0,
            vector_sampling=0,skip_uninfected=False,print_every_n_events=1000):
        """Simulate model using exact Gillespie algorithm.

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
        skip_uninfected -- whether to save only infected hosts/vectors and
            record the number of uninfected host/vectors instead (Boolean,
            default False)
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
            vector_sampling=vector_sampling,
            skip_uninfected=skip_uninfected
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
                dt = self.random.exponential( 1/r_tot ) # time until next event

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
                        dt = self.random.exponential( 1/r_tot )
                            # time until next event
                    else:
                        self.model.t_var = tf # go to end

                self.model.t_var += dt # add time step to main timer

                # Event handling
                if self.model.t_var < tf: # if still within max time
                    u = self.random.random() * r_tot
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
                                                vector_sampling=vector_sampling,
                                                skip_uninfected=skip_uninfected
                                                )
                                            for evt,num in self.model.global_trackers['num_events'].items():
                                                self.model.global_trackers['num_events_over_time'][evt].append(
                                                    self.model.global_trackers['num_events'][evt]
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
            vector_sampling=vector_sampling,
            skip_uninfected=skip_uninfected
            )
        history[tf].history = None
        history[tf].global_trackers = cp.deepcopy( self.model.global_trackers )

        return history


    def runApproximated(self,t0,tf,dt_leap=None,dt_thre=None,time_sampling=0,
            host_sampling=0,vector_sampling=0,skip_uninfected=False,
            print_every_n_events=1000):
        """Simulate model for a specified time using an approximating algorithm.

        Simulates a time series using a mixed Euler-Gillespie algorithm,
        a variation of tau-leaping.

        Arguments:
        t0 -- initial time point to start simulation at (number)
        tf -- initial time point to end simulation at (number)
        dt_leap -- time leap size used to simulate bursts; if None, set to
            minimum growth threshold time across all populations (number,
            default None)
        dt_thre -- time threshold below which bursts are used; if None, set to
            dt_leap (number, default None)
        time_sampling -- how many events to skip before saving a snapshot of the
            system state (saves all by default), if <0, saves only final state
            (int, default 0)
        host_sampling -- how many hosts to skip before saving one in a snapshot
            of the system state (saves all by default) (int, default 0)
        vector_sampling -- how many vectors to skip before saving one in a
            snapshot of the system state (saves all by default) (int, default 0)
        skip_uninfected -- whether to save only infected hosts/vectors and
            record the number of uninfected host/vectors instead (Boolean,
            default False)
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
            vector_sampling=vector_sampling,
            skip_uninfected=skip_uninfected
            ) }
        intervention_tracker = 0
            # keeps track of what the next intervention should be
        self.model.interventions = sorted(
            self.model.interventions, key=lambda i: i.time
            )

        print_counter = 0 # only used to track when to print
        sampling_counter = 0 # used to track when to save a snapshot

        dt_contact_rate_division = 5
        if dt_leap is None:
            dt_leap = min( [
                # min(
                    # p.t_threshold_host, p.t_threshold_vector,
                    1 / max(
                        ( p.contact_rate_host_host +p.contact_rate_host_vector )
                        * dt_contact_rate_division, # arbitrary tolerance
                        1e-6 # arbitrary tolerance
                        )
                    # )
                for p in self.model.populations.values()
                ] )
                    # if leap time step is not specified, set to minimum
                    # intrahost growth threshold time or minimum time before
                    # infection

        dt_leap = dt_leap * 1.0 # make into float

        if dt_thre is None:
            dt_thre = (
                dt_leap / (
                    self.evt_possible.sum() * len( self.model.populations ) *1.0
                         # make into float
                    )
                ) if dt_leap > 0 else 0
                    # if leap threshold is not specified, set to leap length
                    # divided by total number of event types to evaluate
            # dt_thre = dt_leap if dt_thre is None else dt_thre
        # print([ [ #p.t_threshold_host, p.t_threshold_vector,
        #         1 / max(
        #             ( p.contact_rate_host_host +p.contact_rate_host_vector )
        #             * 10, # arbitrary tolerance
        #             1e-6 # arbitrary tolerance
        #             ) ] for p in self.model.populations.values() ])
        # print(dt_leap,dt_thre)
        # if dt_leap == 0:
        #     dt_thre = 0

        while self.model.t_var < tf: # repeat until t reaches end of timecourse
            population_ids = list( self.model.populations.keys() )
            r = self.getRates(population_ids) # get event rates in this state
            r_tot = np.sum(r) # sum of all rates

            # Time handling
            if r_tot > 0:
                dt = self.random.exponential( 1/r_tot ) # time until next event
                if dt < dt_thre: # if time step is smaller than threshold
                    dt = dt_leap # set time step to tau leap size

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

                        dt_leap = min( [
                            # min(
                                # p.t_threshold_host, p.t_threshold_vector,
                                1 / max(
                                    ( p.contact_rate_host_host +p.contact_rate_host_vector )
                                    * dt_contact_rate_division, # arbitrary tolerance
                                    1e-12 # arbitrary tolerance
                                    )
                                # )
                            for p in self.model.populations.values()
                            ] )
                                # if leap time step is not specified, set to minimum
                                # intrahost growth threshold time or minimum time before
                                # infection

                        dt_leap = dt_leap * 1.0 # make into float

                        dt_thre = (
                            dt_leap / (
                                self.evt_possible.sum() * len( self.model.populations ) *1.0
                                     # make into float
                                )
                            ) if dt_leap > 0 else 0

                        # now recalculate rates
                        population_ids = list( self.model.populations.keys() )
                        r = self.getRates(population_ids)
                            # get event rates in this state
                        r_tot = np.sum(r) # sum of all rates

                    if r_tot > 0: # if no more events happening,
                        dt = self.random.exponential( 1/r_tot )
                            # time until next event
                    else:
                        self.model.t_var = tf # go to end

                self.model.t_var += dt # add time step to main timer

                # Event handling
                if self.model.t_var < tf: # if still within max time
                    if dt == dt_leap: # if we are tau leaping,
                        changed = False

                        # print every n events
                        print_counter += 1
                        if print_counter == print_every_n_events:
                            print_counter = 0
                            print(
                                'Simulating time: '
                                + str(self.model.t_var) + ', event: tau leap burst'
                                )
                        for e in range(r.shape[0]): # for every possible event,
                            for p in range(r.shape[1]):
                                # for every possible population,
                                num_evts_burst = int(self.random.poisson( r[e,p] ))
                                    # compute number of events of this type in
                                    # this population

                                for _ in range(num_evts_burst):
                                    changed_tmp = self.doAction(
                                        e, self.model.populations[
                                            population_ids[p]
                                            ], self.random.random()
                                        ) # do corresponding action,
                                          # feed in random number
                                    changed = True if changed else changed_tmp

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
                                        vector_sampling=vector_sampling,
                                        skip_uninfected=skip_uninfected
                                        )
                                    for evt,num in self.model.global_trackers['num_events'].items():
                                        self.model.global_trackers['num_events_over_time'][evt].append(
                                            self.model.global_trackers['num_events'][evt]
                                            )

                    else:
                        u = self.random.random() * r_tot
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
                                                    vector_sampling=vector_sampling,
                                                    skip_uninfected=skip_uninfected
                                                    )
                                                for evt,num in self.model.global_trackers['num_events'].items():
                                                    self.model.global_trackers['num_events_over_time'][evt].append(
                                                        self.model.global_trackers['num_events'][evt]
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
            vector_sampling=vector_sampling,
            skip_uninfected=skip_uninfected
            )
        history[tf].history = None
        history[tf].global_trackers = cp.deepcopy( self.model.global_trackers )

        return history
