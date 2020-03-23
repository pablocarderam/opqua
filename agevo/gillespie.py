#TODO: adapt to this model's events

class StocModel(object):

    """
    Class defines a model's parameters and methods for changing system state
    according to the possible events and simulating a timecourse using the
    Gillespie algorithm.
    """

    def __init__(self):

        '''
        Class constructor defines parameters and declares state variables.
        '''

        super(StocModel, self).__init__() # initialize as parent class object

        # Event IDs
        self.evt_IDs = [ R_F, R_B, E_F, E_B, N_F, N_B, P_F, P_B, S_F, S_B ]
            # event IDs in specific order


    def getRates(self):

        '''
        Calculates event rates according to current system state.
        Returns:
            dictionary with event ID constants as keys and rates as values.
                Includes total rate under 'tot' key.
        '''

        T = self.mu * np.power( np.sin( (self.t_var-self.xi)*np.pi/self.lam ), 2*self.nu )
            # TF concentration at this time

        rates = np.zeros( [ len(self.evt_IDs),N_GEN ] )
            # rate array size of event space

        rates[R_F,:] = (
            ( self.alp * np.power(T,self.h_T) / ( np.power(T,self.h_T) + np.power(self.k_T,self.h_T) ) ) * ( ( self.x_var[E,:] * np.power(self.x_var[N,:],self.h_R) / ( np.power(self.x_var[N,:],self.h_R) + np.power(self.k_R,self.h_R) ) ) + self.kap )
        )

        rates[R_B,:] = (
            self.dlt * self.x_var[R,:]
        )

        rates[E_F,:] = (
            self.eps * self.x_var[N,:] * np.power( (self.d - self.x_var[E,:]),self.h_E ) / ( ( np.power( (self.d - self.x_var[E,:]),self.h_E ) + np.power(self.k_E,self.h_E) ) )
        )

        rates[E_B,:] = (
            self.zet * np.sum(self.x_var[S,:]) * np.power(self.x_var[E,:],self.h_H) / ( ( np.power(self.x_var[E,:],self.h_H) + np.power(self.k_H,self.h_H) ) )
        )

        rates[N_F,:] = (
            self.gam * self.x_var[R,:] * np.power( (self.n_t - np.sum(self.x_var[N,:]) ),self.h_N ) / ( ( np.power( (self.n_t - np.sum(self.x_var[N,:]) ),self.h_N ) + np.power(self.k_N,self.h_N) ) )
        )

        rates[N_B,:] = (
            self.bet * self.x_var[N,:]
        )

        rates[P_F,:] = (
            self.eta * self.x_var[E,:] * np.power(T,self.h_T) * np.power(self.x_var[N,:],self.h_R) / ( ( np.power(T,self.h_T) + np.power(self.k_T,self.h_T) ) * ( np.power(self.x_var[N,:],self.h_R) + np.power(self.k_R,self.h_R) ) )
        )

        rates[P_B,:] = (
            0
        )

        rates[S_F,:] = (
            self.the * (self.d - self.x_var[E,:]) * np.power(T,self.h_T) / ( np.power(T,self.h_T) + np.power(self.k_T,self.h_T) )
        )

        rates[S_B,:] = (
            self.iot * self.x_var[S,:]
        )

        #print(np.sum(rates,1), T)

        # rate_dict = dict(zip( self.evt_IDs, rates )) # save IDs, rates in dict
        # rate_dict['tot'] = rates.sum() # save total sum of rates

        return rates


    def doAction(self,act,gen):

        '''
        Changes system state variables according to act argument passed (must be
        one of the event ID constants)
        Arguments:
            act : int event ID constant - defines action to be taken
        '''

        if act == R_F:
            self.x_var[R,gen] += 1
        elif act == R_B:
            self.x_var[R,gen] -= 1
        elif act == E_F:
            self.x_var[E,gen] += 1
        elif act == E_B:
            self.x_var[E,gen] -= 1
        elif act == N_F:
            self.x_var[N,gen] += 1
            #self.x_var[R,gen] -= 1
        elif act == N_B:
            self.x_var[N,gen] -= 1
            #self.x_var[R,gen] += 1
        elif act == P_F:
            self.x_var[P,gen] += 1
        elif act == P_B:
            self.x_var[P,:] = 0 # RBC lysis!
        elif act == S_F:
            self.x_var[S,gen] += 1
        elif act == S_B:
            self.x_var[S,gen] -= 1


    def gillespie(self):

        '''
        Simulates a time series with time values specified in argument t_vec
        using the Gillespie algorithm. Stops simulation at maximum time or if no
        change in distance has occurred after the time specified in
        max_no_change. Records position values for the given time values in
        x_vec.
        '''

        # Simulation variables
        self.t_var = self.t[0] # keeps track of time
        self.x[:,:,0] = self.x_var # initialize distance at current distance
        i = 0 # keeps track of index within x and t_vec lists
        intervention_tracker = False # keeps track of what the next intervention should be

        while self.t_var < self.t[-1]:
                # repeat until t reaches end of timecourse
            r = self.getRates() # get event rates in this state
            r_tot = np.sum(r) # sum of all rates
            print(self.t_var)
            if not intervention_tracker and self.t_var > 124: # if there are any interventions left and if it is time to make one,
                # self.x_var[R,0] = 0 # get rid of all aslncRNA for gene 0
                # self.x_var[R,1] += 2500 # add aslncRNA for gene 1

                # self.eps = 0 # get rid of all euchromatin markers
                # self.eps = self.eps*10 # add all euchromatin markers
                # self.zet = 0 # get rid of all heterochromatin markers
                # self.dlt = self.dlt/10 # get rid of all aslncRNAses

                intervention_tracker = True # advance the tracker

            if 1/self.max_dt < r_tot: # if average time to next event is less than maximum permitted time step,
                # allow it, calculate probability
                # Time handling
                dt = np.random.exponential( 1/r_tot ) # time until next event
                self.t_var += dt # add time step to main timer
                while i < len(self.t)-1 and self.t[i+1] < self.t_var:
                    # while still within bounds and the new time exceeds the next
                    # time in list,
                    i += 1 # move to next time frame in list
                    self.x[:,:,i] = self.x[:,:,i-1]
                        # fill in the state with the previous frame
                    self.t_l += self.t[i] - self.t[i-1] # advance life cycle timer
                    if self.t_l > self.lam: # if over life cycle time,
                        self.t_l = 0 # restart life cycle timer
                        self.doAction(P_B,-1) # lysis

                # Event handling
                if self.t_var < self.t[-1]: # if still within max time
                    u = np.random.random() * r_tot
                        # random uniform number between 0 (inc) and total rate (exc)
                    r_cum = 0 # cumulative rate
                    for e in range(r.shape[0]): # for every possible event,
                        for g in range(r.shape[1]): # for every possible gene,
                            r_cum += r[e,g] # add this event's rate to cumulative rate
                            if u < r_cum: # if random number is under cumulative rate
                                self.doAction(e,g) # do corresponding action
                                self.x[:,:,i] = self.x_var # record distance in list
                                break # exit event loop


                        else: # if the inner loop wasn't broken,
                            continue # continue outer loop

                        break # otherwise, break outer loop



            else: # if no event happening probabilistically in this max permitted time step,
                self.t_var += self.max_dt # add time step to main timer
                while i < len(self.t)-1 and self.t[i+1] < self.t_var:
                    # while still within bounds and the new time exceeds the next
                    # time in list,
                    i += 1 # move to next time frame in list
                    self.x[:,:,i] = self.x[:,:,i-1]
                        # fill in the distance with the previous frame
                    self.t_l += self.t[i] - self.t[i-1] # advance life cycle timer
                    if self.t_l > self.lam: # if over life cycle time,
                        self.t_l = 0 # restart life cycle timer
                        self.doAction(P_B,-1) # lysis

        self.t = self.t[0:i+1] # keep only time points that have been simulated
        self.x = self.x[:,:,0:i+1] # keep only distances that have been simulated
        #self.x = self.x / self.vol # make concentration
