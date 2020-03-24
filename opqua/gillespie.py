
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
        self.evt_IDs = [ MIGRATE, CONTACT_HOST_HOST, CONTACT_HOST_VECTOR, RECOVER_HOST, RECOVER_VECTOR ]
            # event IDs in specific order

        self.model = model


    def getRates(model):

        '''
        Calculates event rates according to current system state.
        Returns:
            dictionary with event ID constants as keys and rates as values.
                Includes total rate under 'tot' key.
        '''

        rates = np.zeros( [ len(self.evt_IDs), len(model.populations) ] )
            # rate array size of event space

        rates[MIGRATE,:] = \
            np.array( [ len(p.hosts) * p.total_migration_rate for id,p in model.populations.items() ] )

        rates[CONTACT_HOST_HOST,:] = \
            np.array( [ p.host_host_transmission * len(p.hosts) * len(p.hosts) * p.contact_rate_host_host for id,p in model.populations.items() ] )
                # contact rate assumes fixed area--large populations are dense
                # populations, so contact scales linearly with both host and vector
                # populations. If you don't want this to happen, modify the population's
                # contact rate accordingly.

        rates[CONTACT_HOST_VECTOR,:] = \
            np.array( [ p.vector_borne * len(p.hosts) * len(p.vectors) * p.contact_rate_host_vector for id,p in model.populations.items() ] )
                # contact rate assumes fixed area--large populations are dense
                # populations, so contact scales linearly with both host and vector
                # populations. If you don't want this to happen, modify the population's
                # contact rate accordingly.

        rates[RECOVER_HOST,:] = \
            np.array( [ p.num_infected_hosts * p.recovery_rate_host for id,p in model.populations.items() ] )

        rates[RECOVER_VECTOR,:] = \
            np.array( [ p.num_infected_vectors * p.recovery_rate_vector for id,p in model.populations.items() ] )

        #print(np.sum(rates,1), T)

        # rate_dict = dict(zip( self.evt_IDs, rates )) # save IDs, rates in dict
        # rate_dict['tot'] = rates.sum() # save total sum of rates

        return rates


    def doAction(self,act,pop,rand):

        '''
        Changes system state variables according to act argument passed (must be
        one of the event ID constants)
        Arguments:
            act : int event ID constant - defines action to be taken
        '''

        if act == MIGRATE:
            rand = rand * pop.total_migration_rate
            r_cum = 0
            for neighbor in pop.neighbors:
                r_cum += pop.neighbors[neighbor]
                if r_cum > rand:
                    pop.migrate(neighbor,1,0)



        elif act == CONTACT_HOST_VECTOR:
            rand = rand * len(pop.hosts)
            host = np.floor(rand)
            vector = np.floor( ( rand - host ) * len(pop.vectors) )
            pop.contactVectorBorne(host,vector)

        elif act == CONTACT_HOST_HOST:
            rand = rand * len(pop.hosts)
            host1 = np.floor(rand)
            host2 = np.floor( ( rand - host1 ) * len(pop.hosts) )
            pop.contactHostHost(host1,host2)

        elif act == RECOVER_HOST:
            host = np.floor( rand * len(pop.hosts) )
            pop.recoverHost(host)

        elif act == RECOVER_VECTOR:
            vector = np.floor( rand * len(pop.vectors) )
            pop.recoverVector(vector)


    def addToDf(model,df,time):
        """ Saves status of model to dataframe given """

        df = pd.concat([df] +
            [
                [
                    pd.DataFrame( time, pop.id, 'Host', host.id, [ str(host.parasites) ],
                        columns=['Time','Population','Organism','ID','Parasites'] )
                for host in pop.hosts ]
            for id,pop in model.populations.items() ] +
            [
                [
                    pd.DataFrame( time, pop.id, 'Vector', vector.id, [ str(vector.parasites) ],
                        columns=['Time','Population','Organism','ID','Parasites'] )
                for vector in pop.vectors ]
            for id,pop in model.populations.items() ],

          ignore_index=True
        )

        return df


    def run(t0,tf):

        '''
        Simulates a time series with time values specified in argument t_vec
        using the Gillespie algorithm. Stops simulation at maximum time or if no
        change in distance has occurred after the time specified in
        max_no_change. Records position values for the given time values in
        x_vec.
        '''

        # Simulation variables
        t_var = t0 # keeps track of time
        dat = pd.DataFrame( columns=['Time','Population','Organism','ID','Parasites'] )
        intervention_tracker = 0 # keeps track of what the next intervention should be
        model.interventions = sorted(model.interventions, key=lambda i: i.time)

        while t_var < tf:
                # repeat until t reaches end of timecourse
            r = self.getRates() # get event rates in this state
            r_tot = np.sum(r) # sum of all rates

            if intervention_tracker < len(model.interventions) and t_var > model.interventions[intervention_tracker].time: # if there are any interventions left and if it is time to make one,
                model.interventions[intervention_tracker].doIntervention()
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
                            self.doAction( e, p, ( u - r_cum + r[e,p] ) / r[e,p] ) # do corresponding action, feed in renormalized random number
                            dat = self.addToDf(model,dat) # record model state in dataframe
                            break # exit event loop


                    else: # if the inner loop wasn't broken,
                        continue # continue outer loop

                    break # otherwise, break outer loop



        return dat
