
import copy as cp
import pandas as pd
import numpy as np
from opqua.classes import *
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

    def __init__(self, model):

        '''
        Class constructor defines parameters and declares state variables.
        '''

        super(Gillespie, self).__init__() # initialize as parent class object

        # Event IDs
        self.evt_IDs = [ self.MIGRATE, self.CONTACT_HOST_HOST, self.CONTACT_HOST_VECTOR, self.RECOVER_HOST, self.RECOVER_VECTOR ]
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
            np.array( [ self.model.populations[id].host_host_transmission * len(self.model.populations[id].hosts) * len(self.model.populations[id].hosts) * self.model.populations[id].contact_rate_host_host for id,p in self.model.populations.items() ] )
                # contact rate assumes fixed area--large populations are dense
                # populations, so contact scales linearly with both host and vector
                # populations. If you don't want this to happen, modify the population's
                # contact rate accordingly.

        rates[self.CONTACT_HOST_VECTOR,:] = \
            np.array( [ self.model.populations[id].vector_borne * len(self.model.populations[id].hosts) * len(self.model.populations[id].vectors) * self.model.populations[id].contact_rate_host_vector for id,p in self.model.populations.items() ] )
                # contact rate assumes fixed area--large populations are dense
                # populations, so contact scales linearly with both host and vector
                # populations. If you don't want this to happen, modify the population's
                # contact rate accordingly.

        rates[self.RECOVER_HOST,:] = \
            np.array( [ len( self.model.populations[id].infected_hosts ) * self.model.populations[id].recovery_rate_host for id,p in self.model.populations.items() ] )

        rates[self.RECOVER_VECTOR,:] = \
            np.array( [ len( self.model.populations[id].infected_vectors ) * self.model.populations[id].recovery_rate_vector for id,p in self.model.populations.items() ] )

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
            host = int( np.floor( rand * len(pop.hosts) ) )
            pop.recoverHost(host)
            changed = True

        elif act == self.RECOVER_VECTOR:
            vector = int( np.floor( rand * len(pop.vectors) ) )
            pop.recoverVector(vector)
            changed = True

        return changed


    def saveToDf(self,history,df):
        """ Saves status of model to dataframe given """

        print('Saving file...')
        df = pd.concat([df] + #TODO: this can be parallelized with jl
            [ pd.concat(
                [ pd.concat(
                        [ pd.DataFrame( [ [ time, pop.id, 'Host', host.id, str( list( host.pathogens.keys() ) ) ] ],
                            columns=['Time','Population','Organism','ID','Pathogens'] )
                            for host in pop.hosts ] + [
                        pd.DataFrame( [ [ time, pop.id, 'Vector', vector.id, str( list( vector.pathogens.keys() ) ) ] ],
                                    columns=['Time','Population','Organism','ID','Pathogens'] )
                            for vector in pop.vectors ]
                    )
                    for id,pop in model.populations.items() ]
                )
                for time,model in history.items() ],
            ignore_index=True)

        print('...file saved.')

        return df


    def run(self,t0,tf):

        '''
        Simulates a time series with time values specified in argument t_vec
        using the Gillespie algorithm. Stops simulation at maximum time or if no
        change in distance has occurred after the time specified in
        max_no_change. Records position values for the given time values in
        x_vec.
        '''

        # Simulation variables
        t_var = t0 # keeps track of time
        dat = pd.DataFrame( columns=['Time','Population','Organism','ID','Pathogens'] )
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
                            t_inf_hosts = np.array([(len(h.pathogens)>0) for h in self.model.populations['p1'].hosts]).sum()
                            t_inf_vectors = np.array([(len(v.pathogens)>0) for v in self.model.populations['p1'].vectors]).sum()
                            print('Simulating time: ' + str(t_var),e,len(self.model.populations['p1'].infected_hosts), t_inf_hosts, len(self.model.populations['p1'].infected_vectors), t_inf_vectors)
                            changed = self.doAction( e, self.model.populations[ population_ids[p] ], ( u - r_cum + r[e,p] ) / r[e,p] ) # do corresponding action, feed in renormalized random number
                            #dat = self.addToDf(self.model,dat,t_var) # record model state in dataframe
                            if changed:
                                history[t_var] = cp.deepcopy(self.model)

                            break # exit event loop


                    else: # if the inner loop wasn't broken,
                        continue # continue outer loop

                    break # otherwise, break outer loop



        dat = self.saveToDf(history,dat) # record model state in dataframe
        return dat
