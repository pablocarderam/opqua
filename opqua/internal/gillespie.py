# TODO: print death, protection sequences to results
import copy as cp
import pandas as pd
import numpy as np
import joblib as jl

from opqua.internal.classes import *
from opqua.model import *

class Gillespie(object):

    """
    Class defines a model's parameters and methods for changing system state
    according to the possible events and simulating a timecourse using the
    Gillespie algorithm.
    """

    MIGRATE = 0
    CONTACT_HOST_HOST = 1
    CONTACT_HOST_VECTOR = 2
    RECOVER_HOST = 3
    RECOVER_VECTOR = 4
    MUTATE_HOST = 5
    MUTATE_VECTOR = 6
    RECOMBINE_HOST = 7
    RECOMBINE_VECTOR = 8

    def __init__(self, model):

        '''
        Class constructor defines parameters and declares state variables.
        '''

        super(Gillespie, self).__init__() # initialize as parent class object

        # Event IDs
        self.evt_IDs = [ self.MIGRATE, self.CONTACT_HOST_HOST, self.CONTACT_HOST_VECTOR,
            self.RECOVER_HOST, self.RECOVER_VECTOR, self.MUTATE_HOST, self.MUTATE_VECTOR,
            self.RECOMBINE_HOST, self.RECOMBINE_VECTOR ]
            # event IDs in specific order

        self.model = model


    def getRates(self,population_ids):

        '''
        Calculates event rates according to current system state.
        Returns:
            dictionary with event ID constants as keys and rates as values.
                Includes total rate under 'tot' key.
        '''

        rates = np.zeros( [ len(self.evt_IDs), len(population_ids) ] )
            # rate array size of event space

        rates[self.MIGRATE,:] = \
            np.array( [ len(self.model.populations[id].hosts) * self.model.populations[id].total_migration_rate for id in population_ids ] )

        rates[self.CONTACT_HOST_HOST,:] = \
            np.array( [ ( len( self.model.populations[id].infected_hosts ) > 0 ) * len(self.model.populations[id].hosts) * len(self.model.populations[id].hosts) * self.model.populations[id].contact_rate_host_host for id,p in self.model.populations.items() ] )
                # contact rate assumes fixed area--large populations are dense
                # populations, so contact scales linearly with both host and vector
                # populations. If you don't want this to happen, modify the population's
                # contact rate accordingly.

        rates[self.CONTACT_HOST_VECTOR,:] = \
            np.array( [ ( len( self.model.populations[id].infected_hosts ) + len( self.model.populations[id].infected_vectors ) > 0 ) * len(self.model.populations[id].hosts) * len(self.model.populations[id].vectors) * self.model.populations[id].contact_rate_host_vector for id,p in self.model.populations.items() ] )
                # contact rate assumes fixed area--large populations are dense
                # populations, so contact scales linearly with both host and vector
                # populations. If you don't want this to happen, modify the population's
                # contact rate accordingly.

        rates[self.RECOVER_HOST,:] = \
            np.array( [ len( self.model.populations[id].infected_hosts ) * self.model.populations[id].recovery_rate_host for id,p in self.model.populations.items() ] )

        rates[self.RECOVER_VECTOR,:] = \
            np.array( [ len( self.model.populations[id].infected_vectors ) * self.model.populations[id].recovery_rate_vector for id,p in self.model.populations.items() ] )

        rates[self.MUTATE_HOST,:] = \
            np.array( [ len( self.model.populations[id].infected_hosts ) * self.model.populations[id].mutate_in_host for id,p in self.model.populations.items() ] )

        rates[self.MUTATE_VECTOR,:] = \
            np.array( [ len( self.model.populations[id].infected_vectors ) * self.model.populations[id].mutate_in_vector for id,p in self.model.populations.items() ] )

        rates[self.RECOMBINE_HOST,:] = \
            np.array( [ len( self.model.populations[id].infected_hosts ) * self.model.populations[id].recombine_in_host for id,p in self.model.populations.items() ] )

        rates[self.RECOMBINE_VECTOR,:] = \
            np.array( [ len( self.model.populations[id].infected_vectors ) * self.model.populations[id].recombine_in_vector for id,p in self.model.populations.items() ] )

        return rates


    def doAction(self,act,pop,rand):

        '''
        Changes system state variables according to act argument passed (must be
        one of the event ID constants)
        Arguments:
            act : int event ID constant - defines action to be taken
        '''

        changed = False

        if act == self.MIGRATE:
            rand = rand * pop.total_migration_rate
            r_cum = 0
            for neighbor in pop.neighbors:
                r_cum += pop.neighbors[neighbor]
                if r_cum > rand:
                    pop.migrate(neighbor,1,0)
                    changed = True



        elif act == self.CONTACT_HOST_VECTOR:
            rand = rand * len(pop.hosts)
            host = int( np.floor(rand) )
            vector = int( np.floor( ( rand - host ) * len(pop.vectors) ) )
            changed = pop.contactVectorBorne(host,vector)

        elif act == self.CONTACT_HOST_HOST:
            rand = rand * len(pop.hosts)
            host1 = int( np.floor(rand) )
            host2 = int( np.floor( ( rand - host1 ) * len(pop.hosts) ) )
            changed = pop.contactHostHost(host1,host2)

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
            pop.mutate(host)
            changed = True

        elif act == self.MUTATE_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.mutate(vector)
            changed = True

        elif act == self.RECOMBINE_HOST:
            host = int( np.floor( rand * len(pop.infected_hosts) ) )
            pop.recombine(host)
            changed = True

        elif act == self.RECOMBINE_VECTOR:
            vector = int( np.floor( rand * len(pop.infected_vectors) ) )
            pop.recombine(vector)
            changed = True

        return changed


    def saveToDf(self,history,save_to_file,n_cores=0):
        """ Saves status of model to dataframe given """

        print('Saving file...')

        if not n_cores:
            n_cores = jl.cpu_count()

        new_df = ','.join( ['Time','Population','Organism','ID','Pathogens'] ) + '\n' + \
            '\n'.join(
                jl.Parallel(n_jobs=n_cores, verbose=10) (
                    jl.delayed( lambda d: ''.join(d) )
                        (
                            '\n'.join( [
                                    '\n'.join( [ ','.join( [ str(time), str(pop.id), 'Host', str(host.id), ';'.join( host.pathogens.keys() ) ] )
                                        for host in pop.hosts ] ) + '\n' +
                                    '\n'.join( [ ','.join( [ str(time), str(pop.id), 'Vector', str(vector.id), ';'.join( vector.pathogens.keys() ) ] )
                                        for vector in pop.vectors ] )
                                for id,pop in model.populations.items()
                            ] )
                        )
                    for time,model in history.items()
                )
            )

        file = open(save_to_file,'w')
        file.write(new_df)
        file.close()

        new_df = pd.read_csv(save_to_file)

        print('...file saved.')

        return new_df


    def run(self,t0,tf,save_to_file):

        '''
        Simulates a time series with time values specified in argument t_vec
        using the Gillespie algorithm. Stops simulation at maximum time or if no
        change in distance has occurred after the time specified in
        max_no_change. Records position values for the given time values in
        x_vec.
        '''

        # Simulation variables
        t_var = t0 # keeps track of time
        history = { 0: cp.deepcopy(self.model) }
        intervention_tracker = 0 # keeps track of what the next intervention should be
        self.model.interventions = sorted(self.model.interventions, key=lambda i: i.time)

        while t_var < tf:
                # repeat until t reaches end of timecourse
            population_ids = list( self.model.populations.keys() )
            r = self.getRates(population_ids) # get event rates in this state
            r_tot = np.sum(r) # sum of all rates

            if intervention_tracker < len(self.model.interventions) and t_var > self.model.interventions[intervention_tracker].time: # if there are any interventions left and if it is time to make one,
                self.model.interventions[intervention_tracker].doIntervention()
                intervention_tracker += 1 # advance the tracker

            # Time handling
            if r_tot > 0:
                dt = np.random.exponential( 1/r_tot ) # time until next event
                t_var += dt # add time step to main timer

                # Event handling
                if t_var < tf: # if still within max time
                    u = np.random.random() * r_tot
                        # random uniform number between 0 (inc) and total rate (exc)
                    r_cum = 0 # cumulative rate
                    for e in range(r.shape[0]): # for every possible event,
                        for p in range(r.shape[1]): # for every possible population,
                            r_cum += r[e,p] # add this event's rate to cumulative rate
                            if u < r_cum: # if random number is under cumulative rate
                                print( 'Simulating time: ' + str(t_var), e, len(self.model.populations['my_population'].infected_hosts), len(self.model.populations['my_population'].infected_vectors) )
                                changed = self.doAction( e, self.model.populations[ population_ids[p] ], ( u - r_cum + r[e,p] ) / r[e,p] ) # do corresponding action, feed in renormalized random number
                                if changed:
                                    history[t_var] = cp.deepcopy(self.model)

                                break # exit event loop


                        else: # if the inner loop wasn't broken,
                            continue # continue outer loop

                        break # otherwise, break outer loop



            else:
                t_var = tf # TODO: next intervention



        history[tf] = cp.deepcopy(self.model)
        dat = self.saveToDf(history,save_to_file) # record model state in dataframe

        return dat
