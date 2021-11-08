
"""Contains class Population."""

import numpy as np
import copy as cp
import random

from opqua.internal.host import Host
from opqua.internal.vector import Vector

class Population(object):
    """Class defines a population with hosts, vectors, and specific parameters.

    Constants:
    # These all denote positions in coefficients_hosts and coefficients_vectors
    INFECTED -- position of "infected" Boolean values for each individual inside
        coefficients array
    CONTACT -- position of intra-population aggregated contact rate for each
        individual inside coefficients array
    LETHALITY -- position of aggregated death rate for each individual inside
        coefficients array
    NATALITY -- position of aggregated birth rate for each individual inside
        coefficients array
    RECOVERY -- position of aggregated recovery rate for each individual inside
        coefficients array
    MIGRATION -- position of aggregated inter-population migration rate for each
        individual inside coefficients array
    POPULATION_CONTACT -- position of inter-population aggregated contact rate
        for each individual inside coefficients array
    MUTATION -- position of aggregated mutation rate for each individual inside
        coefficients array
    RECOMBINATION -- position of aggregated recovery rate for each individual
        inside coefficients array
    NUM_COEFFICIENTS -- total number of types of coefficients (columns) in
        coefficient arrays
    CHROMOSOME_SEPARATOR -- character reserved to denote separate chromosomes in
        genomes

    Methods:
    copyState -- returns a slimmed-down version of the current population state
    setSetup -- assigns a given set of parameters to this population
    addHosts -- adds hosts to the population
    addVectors -- adds vectors to the population
    newHostGroup -- returns a list of random (healthy or any) hosts
    newVectorGroup -- returns a list of random (healthy or any) vectors
    removeHosts -- removes hosts from the population
    removeVectors -- removes vectors from the population
    addPathogensToHosts -- adds pathogens with specified genomes to hosts
    addPathogensToVectors -- adds pathogens with specified genomes to vectors
    treatHosts -- removes infections susceptible to given treatment from hosts
    treatVectors -- removes infections susceptible to treatment from vectors
    protectHosts -- adds protection sequence to hosts
    protectVectors -- adds protection sequence to vectors
    wipeProtectionHosts -- removes all protection sequences from hosts
    wipeProtectionVectors -- removes all protection sequences from hosts
    setHostMigrationNeighbor -- sets migration rate of hosts from this
        population towards another
    setVectorMigrationNeighbor -- sets migration rate of vectors from this
        population towards another
    migrate -- transfers hosts and/or vectors from this population to a neighbor
    setHostPopulationContactNeighbor -- set host contact rate from this
        population towards another one
    setVectorPopulationContactNeighbor -- set vector contact rate from this
        population towards another one
    populationContact -- contacts hosts and/or vectors from this population to
        another
    contactHostHost -- contact any two (weighted) random hosts in population
    contactHostVector -- contact a (weighted) random host and vector in
        population
    contactVectorHost -- contact a (weighted) random vector and host in
        population
    recoverHost --removes all infections from given host
    recoverVector -- removes all infections from given vector
    killHost -- add host at this index to dead list, remove it from alive ones
    killVector -- add vector at this index to dead list, remove it from alive
        ones
    dieHost -- remove host at this index from alive lists
    dieVector -- remove vector at this index from alive lists
    birthHost -- add host at this index to population, remove it from alive ones
    birthVector -- add vector at this index to population, remove it from alive
        ones
    mutateHost -- mutates a single locus in a random pathogen in a host
    mutateVector -- mutates a single locus in a random pathogen in a vector
    recombineHost -- recombines two random pathogens in a host
    recombineVector -- recombines two random pathogens in a host
    updateHostCoefficients -- updates event coefficients in population's hosts
    updateVectorCoefficients -- updates event coefficients in population's
        vectors
    getWeightedRandom -- returns index of element chosen using weights from
        coefficient array and a given random number
    """

    INFECTED = 0
    CONTACT = 1
    RECEIVE_CONTACT = 2
    LETHALITY = 3
    NATALITY = 4
    RECOVERY = 5
    MIGRATION = 6
    POPULATION_CONTACT = 7
    RECEIVE_POPULATION_CONTACT = 8
    MUTATION = 9
    RECOMBINATION = 10

    NUM_COEFFICIENTS = 11

    CHROMOSOME_SEPARATOR = '/'

    def __init__(self, model, id, setup, num_hosts, num_vectors, slim=False):
        """Create a new Population.

        Arguments:
        model -- parent model this population is a part of (Model)
        id -- unique identifier for this population in the model (String)
        setup -- setup object with parameters for this population (Setup)
        num_hosts -- number of hosts to initialize population with (int)
        num_vectors -- number of hosts to initialize population with (int)

        Keyword arguments:
        slim -- whether to create a slimmed-down representation of the
            population for data storage (only ID, host and vector lists)
            (Boolean, default False)
        """
        super(Population, self).__init__()

        self.model = model
        self.id = id

        if not slim:
                # if not slimmed down for data storage, save other attributes

                # Each element in these following arrays contains the sum of all
                # pathogen rate coefficients within a host, weighted by fitness
                # and normalized to sum_fitness to obtain coefficient that
                # modifies the respective population-wide rate. Order given by
                # constants in Population class. Each host is one row.
            self.coefficients_hosts = np.zeros((1,self.NUM_COEFFICIENTS))
                # all weighted event rate coefficient modifiers for hosts
                # in population, first row is a dummy placeholder row
            self.coefficients_vectors = np.zeros((1,self.NUM_COEFFICIENTS))
                # all weighted event rate coefficient modifiers for vectors
                # in population, first row is a dummy placeholder row

            self.hosts = [
                Host(
                    self, self.id + '_' + str(id)
                ) for id in range(int(num_hosts))
                ]
                # contains all live hosts
            self.vectors = [
                Vector(
                    self, self.id + '_' + str(id)
                ) for id in range(int(num_vectors))
                ]
                # contains all live vectors

            # delete dummy first row, correct indexes of hosts afterwards
            # self.coefficients_hosts = np.delete(
            #     self.coefficients_hosts, 0, 0
            #     ) if num_hosts > 0 else self.coefficients_hosts
            # self.coefficients_vectors = np.delete(
            #     self.coefficients_vectors, 0, 0
            #     ) if num_vectors > 0 else self.coefficients_vectors
            # for h in self.hosts:
            #     h.coefficient_index -= 1
            # for v in self.vectors:
            #     v.coefficient_index -= 1

            self.infected_hosts = []
                # contains live, infected hosts (also in hosts list)
            self.healthy_hosts = self.hosts.copy()
                # contains live, healthy hosts (also in hosts list)
            self.infected_vectors = []
                # contains live, infected vectors (also in vectors list)
            self.dead_hosts = [] # contains dead hosts (not in hosts list)
            self.dead_vectors = [] # contains dead vectors (not in vectors list)

            self.neighbors_hosts = {} # dictionary with neighboring populations,
                # keys=Population, values=migration rate from this population to
                # neighboring population, for hosts only
            self.neighbors_vectors = {} # dictionary with neighbor populations,
                # keys=Population, values=migration rate from this population to
                # neighboring population, for vectors only
            self.total_migration_rate_hosts = 0 # sum of all migration rates
                # from this population to neighbors, hosts only
            self.total_migration_rate_vectors = 0 # sum of all migration rates
                # from this population to neighbors, vectors only

            self.neighbors_contact_hosts = {} # dictionary with neighboring
                # populations, keys=Population, values=population contact rate
                # from this population to neighboring population, for hosts only
            self.neighbors_contact_vectors = {} # dictionary with neighbor
                # populations, keys=Population, values=population contact rate
                # from this population to neighboring population, for vectors only
            self.population_contact_rate_host = 0 # weighted sum of all host
                # population contact rates from this population to neighbors
            self.population_contact_rate_vector = 0 # weighted sum of all vector
                # population contact rates from this population to neighbors

            self.setSetup(setup)

            self.setHostMigrationNeighbor(self,0)
            self.setVectorMigrationNeighbor(self,0)
            self.setHostPopulationContactNeighbor(self,0)
            self.setVectorPopulationContactNeighbor(self,0)

    def copyState(self,host_sampling=0,vector_sampling=0):
        """Returns a slimmed-down version of the current population state.

        Arguments:
        host_sampling -- how many hosts to skip before saving one in a snapshot
            of the system state (saves all by default) (int, default 0)
        vector_sampling -- how many vectors to skip before saving one in a
            snapshot of the system state (saves all by default) (int, default 0)

        Returns:
        Population object with current host and vector lists.
        """

        copy = Population(self.model, self.id, None, 0, 0, slim=True)
        if host_sampling > 0:
            host_sample = random.sample(
                self.hosts, int( len(self.hosts) / (host_sampling+1) )
                )
            dead_host_sample = random.sample(
                self.dead_hosts, int( len(self.dead_hosts) / (host_sampling+1) )
                )
            copy.hosts = [ h.copyState() for h in host_sample ]
            copy.dead_hosts = [ h.copyState() for h in dead_host_sample ]
        else:
            copy.hosts = [ h.copyState() for h in self.hosts ]
            copy.dead_hosts = [ h.copyState() for h in self.dead_hosts ]

        if vector_sampling > 0:
            vector_sample = random.sample(
                self.vectors, int( len(self.vectors) / (vector_sampling+1) )
                )
            dead_vector_sample = random.sample(
                self.dead_vectors,int(len(self.dead_vectors)/(vector_sampling+1))
                )
            copy.vectors = [ v.copyState() for v in vector_sample ]
            copy.dead_vectors = [ v.copyState() for v in dead_vector_sample ]
        else:
            copy.vectors = [ v.copyState() for v in self.vectors ]
            copy.dead_vectors = [ v.copyState() for v in self.dead_vectors ]

        return copy

    def setSetup(self, setup):
        """Assign parameters stored in Setup object to this population.

        Arguments:
        setup -- the setup to be assigned (Setup)
        """

        self.setup = setup

        self.num_loci = setup.num_loci
        self.possible_alleles = setup.possible_alleles

        self.fitnessHost = setup.fitnessHost
        self.contactHost = setup.contactHost
        self.receiveContactHost = setup.receiveContactHost
        self.lethalityHost = setup.lethalityHost
        self.natalityHost = setup.natalityHost
        self.recoveryHost = setup.recoveryHost
        self.migrationHost = setup.migrationHost
        self.populationContactHost = setup.populationContactHost
        self.receivePopulationContactHost = setup.receivePopulationContactHost
        self.mutationHost = setup.mutationHost
        self.recombinationHost = setup.recombinationHost
        self.updateHostCoefficients()

        self.fitnessVector = setup.fitnessVector
        self.contactVector = setup.contactVector
        self.receiveContactVector = setup.receiveContactVector
        self.lethalityVector = setup.lethalityVector
        self.natalityVector = setup.natalityVector
        self.recoveryVector = setup.recoveryVector
        self.migrationVector = setup.migrationVector
        self.populationContactVector = setup.populationContactVector
        self.receivePopulationContactVector = setup.receivePopulationContactVector
        self.mutationVector = setup.mutationVector
        self.recombinationVector = setup.recombinationVector
        self.updateVectorCoefficients()

        self.contact_rate_host_vector = setup.contact_rate_host_vector
        self.contact_rate_host_host = setup.contact_rate_host_host
            # contact rate assumes fixed area--large populations are dense
            # populations, so contact scales linearly with both host and vector
            # populations. If you don't want this to happen, modify the
            # population's contact rate accordingly.
        self.transmission_efficiency_host_vector = setup.transmission_efficiency_host_vector
        self.transmission_efficiency_vector_host = setup.transmission_efficiency_vector_host
        self.transmission_efficiency_host_host = setup.transmission_efficiency_host_host
        self.mean_inoculum_host = setup.mean_inoculum_host
        self.mean_inoculum_vector = setup.mean_inoculum_vector
        self.recovery_rate_host = setup.recovery_rate_host
        self.recovery_rate_vector = setup.recovery_rate_vector
        self.lethality_rate_host = setup.lethality_rate_host
        self.lethality_rate_vector = setup.lethality_rate_vector
        self.recombine_in_host = setup.recombine_in_host
        self.recombine_in_vector = setup.recombine_in_vector
        self.num_crossover_host = setup.num_crossover_host
        self.num_crossover_vector = setup.num_crossover_vector
        self.mutate_in_host = setup.mutate_in_host
        self.mutate_in_vector = setup.mutate_in_vector
        self.death_rate_host = setup.death_rate_host
        self.death_rate_vector = setup.death_rate_vector
        self.birth_rate_host = setup.birth_rate_host
        self.birth_rate_vector = setup.birth_rate_vector
        self.vertical_transmission_host = setup.vertical_transmission_host
        self.vertical_transmission_vector = setup.vertical_transmission_vector
        self.inherit_protection_host = setup.inherit_protection_host
        self.inherit_protection_vector = setup.inherit_protection_vector
        self.protection_upon_recovery_host \
            = setup.protection_upon_recovery_host
        self.protection_upon_recovery_vector \
            = setup.protection_upon_recovery_vector

    def addHosts(self, num_hosts):
        """Add a number of healthy hosts to population, return list with them.

        Arguments:
        num_hosts -- number of hosts to be added (int)

        Returns:
        list containing new hosts
        """

        new_hosts = [
            Host(
                self, self.id + '_' + str( i + len(self.hosts) )
                ) for i in range(num_hosts)
            ]
        self.hosts += new_hosts
        self.healthy_hosts += new_hosts

        return new_hosts

    def addVectors(self, num_vectors):
        """Add a number of healthy vectors to population, return list with them.

        Arguments:
        num_vectors -- number of vectors to be added (int)

        Returns:
        list containing new vectors
        """

        new_vectors = [
            Vector(
                self, self.id + '_' + str( i + len(self.vectors) )
                ) for i in range(num_vectors)
            ]
        self.vectors += new_vectors

        return new_vectors

    def newHostGroup(self, hosts=-1, type='any'):
        """Return a list of random hosts in population.

        Arguments:
        hosts -- number of hosts to be sampled randomly: if <0, samples from
            whole population; if <1, takes that fraction of population; if >=1,
            samples that integer number of hosts (default -1, number)

        Keyword arguments:
        type -- whether to sample healthy hosts only, infected hosts only, or
            any hosts (default 'any'; String = {'healthy', 'infected', 'any'})

        Returns:
        list containing sampled hosts
        """

        possible_hosts = []

        if type=='healthy':
            possible_hosts = self.healthy_hosts
        elif type=='infected':
            possible_hosts = self.infected_hosts
        elif type=='any':
            possible_hosts = self.hosts
        else:
            raise ValueError(
                '"' + str(type)
                + '" is not a type of host for newHostGroup.'
                )

        num_hosts = -1
        if hosts < 0:
            num_hosts = len(possible_hosts)
        elif hosts < 1:
            num_hosts = int( hosts * len(possible_hosts) )
        else:
            num_hosts = hosts

        if len(possible_hosts) >= num_hosts:
            return np.random.choice(possible_hosts, num_hosts, replace=False)
        else:
            raise ValueError(
                "You're asking for " + str(num_hosts)
                + '"' + type + '" hosts, but population ' + str(self.id)
                + " only has " + str( len(self.healthy_hosts) ) + "."
                )

    def newVectorGroup(self, vectors=-1, type='any'):
        """Return a list of random vectors in population.

        Arguments:
        vectors -- number of vectors to be sampled randomly: if <0, samples from
            whole population; if <1, takes that fraction of population; if >=1,
            samples that integer number of vectors (default -1, number)

        Keyword arguments:
        type -- whether to sample healthy vectors only, infected vectors
            only, or any vectors (default 'any'; String = {'healthy',
            'infected', 'any'})

        Returns:
        list containing sampled vectors
        """

        possible_vectors = []
        if type=='healthy':
            for vector in self.vectors:
                if vector not in self.infected_vectors:
                    possible_vectors.append(vector)
        elif type=='infected':
            possible_vectors = self.infected_vectors
        elif type=='any':
            possible_vectors = self.vectors
        else:
            raise ValueError(
                '"' + str(type)
                + '" is not a type of vector for newVectorGroup.'
                )

        num_vectors = -1
        if vectors < 0:
            num_vectors = len(possible_vectors)
        elif vectors < 1:
            num_vectors = int( vectors * len(possible_vectors) )
        else:
            num_vectors = vectors

        if len(possible_vectors) >= num_vectors:
            return np.random.choice(possible_vectors, num_vectors, replace=False)
        else:
            raise ValueError(
                "You're asking for " + str(num_vectors)
                + '"' + type + '" vectors, but population ' + str(self.id)
                + " only has "
                + str( len(self.vectors) - len(self.infected_vectors) )
                + "."
                )

    def removeHosts(self, num_hosts_or_list):
        """Remove a number of specified or random hosts from population.

        Arguments:
        num_hosts_or_list -- number of hosts to be sampled randomly for removal
            or list of hosts to be removed, must be hosts in this population
            (int or list of Hosts)
        """

        if isinstance(num_hosts_or_list, list):
            for host_removed in num_hosts_or_list:
                if host_removed in self.hosts:
                    if host_removed in self.infected_hosts:
                        self.infected_hosts.remove( host_removed )
                    else:
                        self.healthy_hosts.remove( host_removed )

                    self.hosts.remove( host_removed )
                    for h in self.hosts:
                        if h.coefficient_index > host_removed.coefficient_index:
                            h.coefficient_index -= 1

                    self.coefficients_hosts = np.delete(
                        self.coefficients_hosts,
                        host_removed.coefficient_index, 0
                        )
        else:
            for _ in range(num_hosts_or_list):
                host_removed = np.random.choice(self.hosts)
                if host_removed in self.infected_hosts:
                    self.infected_hosts.remove( host_removed )
                else:
                    self.healthy_hosts.remove( host_removed )

                self.hosts.remove( host_removed )
                for h in self.hosts:
                    if h.coefficient_index > host_removed.coefficient_index:
                        h.coefficient_index -= 1

                self.coefficients_hosts = np.delete(
                    self.coefficients_hosts,
                    host_removed.coefficient_index, 0
                    )

    def removeVectors(self, num_vectors_or_list):
        """Remove a number of specified or random vectors from population.

        Arguments:
        num_vectors_or_list -- number of vectors to be sampled randomly for
            removal or list of vectors to be removed, must be vectors in this
            population (int or list of Vectors)
        """

        if isinstance(num_vectors_or_list, list):
            for vector_removed in num_vectors_or_list:
                if vector_removed in self.vectors:
                    if vector_removed in self.infected_vectors:
                        self.infected_vectors.remove( vector_removed )

                    self.vectors.remove( vector_removed )
                    for v in self.vectors:
                        if v.coefficient_index>vector_removed.coefficient_index:
                            v.coefficient_index -= 1

                    self.coefficients_vectors = np.delete(
                        self.coefficients_vectors,
                        vector_removed.coefficient_index, 0
                        )
        else:
            for _ in range(num_vectors):
                vector_removed = np.random.choice(self.vectors)
                if vector_removed in self.infected_vectors:
                    self.infected_vectors.remove( vector_removed )

                self.vectors.remove( vector_removed )
                for v in self.vectors:
                    if v.coefficient_index > vector_removed.coefficient_index:
                        v.coefficient_index -= 1

                self.coefficients_vectors = np.delete(
                    self.coefficients_vectors,
                    vector_removed.coefficient_index, 0
                    )

    def addPathogensToHosts(self, genomes_numbers, hosts=[]):
        """Add specified pathogens to random hosts, optionally from a list.

        Arguments:
        genomes_numbers -- dictionary conatining pathogen genomes to add as keys
            and number of hosts each one will be added to as values (dict with
            keys=Strings, values=int)

        Keyword arguments:
        hosts -- list of specific hosts to sample from, if empty, samples from
            whole population (default empty list; empty)
        """

        if len(hosts) == 0:
            hosts = self.hosts

        for genome in genomes_numbers:
            if len(genome) == self.num_loci and genomes_numbers[genome]>0 and all( [
                allele in self.possible_alleles[i] + self.CHROMOSOME_SEPARATOR
                    for i,allele in enumerate(genome)
                ] ):
                rand_hosts = np.random.choice(
                    hosts, int(genomes_numbers[genome]), replace=False
                    )

                for host in rand_hosts:
                    host.acquirePathogen(genome)

                if genome not in self.model.global_trackers['genomes_seen']:
                    self.model.global_trackers['genomes_seen'].append(genome)

            elif genomes_numbers[genome]>0:
                raise ValueError('Genome ' + genome + ' must be of length '
                    + str(self.num_loci)
                    + ' and contain only the following characters at each '
                    + 'position: ' + str(self.possible_alleles)
                    + ', or the chromosome separator character "'
                    + self.CHROMOSOME_SEPARATOR + '" .')

    def addPathogensToVectors(self, genomes_numbers, vectors=[]):
        """Add specified pathogens to random vectors, optionally from a list.

        Arguments:
        genomes_numbers -- dictionary conatining pathogen genomes to add as keys
            and number of vectors each one will be added to as values (dict with
            keys=Strings, values=int)

        Keyword arguments:
        vectors -- list of specific vectors to sample from, if empty, samples
            from whole population (default empty list; empty)
        """

        if len(vectors) == 0:
            vectors = self.vectors

        for genome in genomes_numbers:
            if len(genome) == self.num_loci and genomes_numbers[genome]>0 and all( [
                allele in self.possible_alleles[i] + self.CHROMOSOME_SEPARATOR
                    for i,allele in enumerate(genome)
                ] ):
                rand_vectors = np.random.choice(
                    vectors, int(genomes_numbers[genome]), replace=False
                    )

                for vector in rand_vectors:
                    vector.acquirePathogen(genome)

                if genome not in self.model.global_trackers['genomes_seen']:
                    self.model.global_trackers['genomes_seen'].append(genome)

            elif genomes_numbers[genome]>0:
                raise ValueError('Genome ' + genome + ' must be of length '
                    + str(self.num_loci)
                    + ' and contain only the following characters at each '
                    + 'position: ' + self.possible_alleles
                    + self.CHROMOSOME_SEPARATOR + ' .')

    def treatHosts(self, frac_hosts, resistance_seqs, hosts=[]):
        """Treat random fraction of infected hosts against some infection.

        Removes all infections with genotypes susceptible to given treatment.
        Pathogens are removed if they are missing at least one of the sequences
        in resistance_seqs from their genome. Removes this organism from
        population infected list and adds to healthy list if appropriate.

        Arguments:
        frac_hosts -- fraction of hosts considered to be randomly selected
            (number between 0 and 1)
        resistance_seqs -- contains sequences required for treatment resistance
            (list of Strings)

        Keyword arguments:
        hosts -- list of specific hosts to sample from, if empty, samples from
            whole population (default empty list; empty)
        """

        hosts_to_consider = self.hosts
        if len(hosts) > 0:
            hosts_to_consider = hosts

        possible_infected_hosts = []
        for host in hosts_to_consider:
            if len( host.pathogens ):
                possible_infected_hosts.append( host )

        treat_hosts = np.random.choice(
            possible_infected_hosts,
            int( frac_hosts * len( possible_infected_hosts ) ), replace=False
            )
        for host in treat_hosts:
            host.applyTreatment(resistance_seqs)

    def treatVectors(self, frac_vectors, resistance_seqs, vectors=[]):
        """Treat random fraction of infected vectors agains some infection.

        Removes all infections with genotypes susceptible to given treatment.
        Pathogens are removed if they are missing at least one of the sequences
        in resistance_seqs from their genome. Removes this organism from
        population infected list and adds to healthy list if appropriate.

        Arguments:
        frac_vectors -- fraction of vectors considered to be randomly selected
            (number between 0 and 1)
        resistance_seqs -- contains sequences required for treatment resistance
            (list of Strings)

        Keyword arguments:
        vectors -- list of specific vectors to sample from, if empty, samples
            from whole population (default empty list; empty)
        """

        vectors_to_consider = self.vectors
        if len(vectors) > 0:
            vectors_to_consider = vectors

        possible_infected_vectors = []
        for vector in vectors_to_consider:
            if len( vector.pathogens ):
                possible_infected_vectors.append( vector )

        treat_vectors = np.random.choice(
            possible_infected_vectors,
            int( frac_vectors * len( possible_infected_vectors ) ),
            replace=False
            )
        for vector in treat_vectors:
            vector.applyTreatment(resistance_seqs)

    def protectHosts(self, frac_hosts, protection_sequence, hosts=[]):
        """Protect a random fraction of infected hosts against some infection.

        Adds protection sequence specified to a random fraction of the hosts
        specified. Does not cure them if they are already infected.

        Arguments:
        frac_hosts -- fraction of hosts considered to be randomly selected
            (number between 0 and 1)
        protection_sequence -- sequence against which to protect (String)

        Keyword arguments:
        hosts -- list of specific hosts to sample from, if empty, samples from
            whole population (default empty list; empty)
        """

        hosts_to_consider = self.hosts
        if len(hosts) > 0:
            hosts_to_consider = hosts

        protect_hosts = np.random.choice(
            self.hosts, int( frac_hosts * len( hosts_to_consider ) ),
            replace=False
            )
        for host in protect_hosts:
            host.protection_sequences.append(protection_sequence)

    def protectVectors(self, frac_vectors, protection_sequence, vectors=[]):
        """Protect a random fraction of infected vectors against some infection.

        Adds protection sequence specified to a random fraction of the vectors
        specified. Does not cure them if they are already infected.

        Arguments:
        frac_vectors -- fraction of vectors considered to be randomly selected
            (number between 0 and 1)
        protection_sequence -- sequence against which to protect (String)

        Keyword arguments:
        vectors -- list of specific vectors to sample from, if empty, samples
            from whole population (default empty list; empty)
        """

        vectors_to_consider = self.vectors
        if len(vectors) > 0:
            vectors_to_consider = vectors

        protect_vectors = np.random.choice(
            self.vectors, int( frac_vectors * len( vectors_to_consider ) ),
            replace=False
            )
        for vector in protect_vectors:
            vector.protection_sequences.append(protection_sequence)

    def wipeProtectionHosts(self, hosts=[]):
        """Removes all protection sequences from hosts.

        Keyword arguments:
        hosts -- list of specific hosts to sample from, if empty, samples from
            whole population (default empty list; empty)
        """

        hosts_to_consider = self.hosts
        if len(hosts) > 0:
            hosts_to_consider = hosts

        for host in hosts_to_consider:
            host.protection_sequences = []

    def wipeProtectionVectors(self, vectors=[]):
        """Removes all protection sequences from vectors.

        Keyword arguments:
        vectors -- list of specific vectors to sample from, if empty, samples from
            whole population (default empty list; empty)
        """

        vectors_to_consider = self.vectors
        if len(vectors) > 0:
            vectors_to_consider = vectors

        for vector in vectors_to_consider:
            vector.protection_sequences = []

    def setHostMigrationNeighbor(self, neighbor, rate):
        """Set host migration rate from this population towards another one.

         Arguments:
         neighbor -- population towards which migration rate will be specified
            (Population)
         rate -- migration rate from this population to the neighbor (number)
         """

        if neighbor in self.neighbors_hosts:
            self.total_migration_rate_hosts -= self.neighbors_hosts[neighbor]

        self.neighbors_hosts[neighbor] = rate
        self.total_migration_rate_hosts += rate

    def setVectorMigrationNeighbor(self, neighbor, rate):
        """Set vector migration rate from this population towards another one.

         Arguments:
         neighbor -- population towards which migration rate will be specified
            (Population)
         rate -- migration rate from this population to the neighbor (number)
         """

        if neighbor in self.neighbors_vectors:
            self.total_migration_rate_vectors -= self.neighbors_vectors[neighbor]

        self.neighbors_vectors[neighbor] = rate
        self.total_migration_rate_vectors += rate

    def migrate(self, target_pop, num_hosts, num_vectors, rand=None):
        """Transfer hosts and/or vectors from this population to another.

        Arguments:
        target_pop -- population towards which migration will occur (Population)
        num_hosts -- number of hosts to transfer (int)
        num_vectors -- number of vectors to transfer (int)

        Keyword arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individuals to migrate; if None, generates new random number to
            choose (through numpy), otherwise, assumes event is happening
            through Gillespie class call and migrates a single host or vector
            (default None; 0-1 number)
        """

        if rand is None:
            migrating_hosts = np.random.choice(
                self.hosts, num_hosts, replace=False,
                p=self.coefficients_hosts[:,self.MIGRATION]
                )
            migrating_vectors = np.random.choice(
                self.vectors, num_vectors, replace=False,
                p=self.coefficients_vectors[:,self.MIGRATION]
                )
        elif num_hosts == 1:
            index_host,rand = self.getWeightedRandom(
                rand, self.coefficients_hosts[:,self.MIGRATION]
                )
            migrating_hosts = [self.hosts[index_host]]
            migrating_vectors = []
        else:
            index_vector,rand = self.getWeightedRandom(
                rand, self.coefficients_vectors[:,self.MIGRATION]
                )
            migrating_hosts = []
            migrating_vectors = [self.vectors[index_vector]]


        for host in migrating_hosts:
            genomes = { g:1 for g in host.pathogens }
            self.removeHosts([host])
            new_host_list = target_pop.addHosts(1) # list of 1
            target_pop.addPathogensToHosts(genomes,hosts=new_host_list)
            host = None

        for vector in migrating_vectors:
            genomes = { g:1 for g in vector.pathogens }
            self.removeVectors([vector])
            new_vector_list = target_pop.addVectors(1) # list of 1
            target_pop.addPathogensToVectors(genomes,vectors=new_vector_list)
            vector = None

    def setHostPopulationContactNeighbor(self, neighbor, rate):
        """Set host contact rate from this population towards another one.

         Arguments:
         neighbor -- population towards which migration rate will be specified
            (Population)
         rate -- migration rate from this population to the neighbor (number)
         """

        self.neighbors_contact_hosts[neighbor] = rate

    def setVectorPopulationContactNeighbor(self, neighbor, rate):
        """Set vector contact rate from this population towards another one.

         Arguments:
         neighbor -- population towards which migration rate will be specified
            (Population)
         rate -- migration rate from this population to the neighbor (number)
         """

        self.neighbors_contact_vectors[neighbor] = rate

    def populationContact(
        self, target_pop, rand, host_origin=True, host_target=True):
        """Contacts hosts and/or vectors from this population to another.

        Arguments:
        target_pop -- population towards which migration will occur (Population)
        rand -- uniform random number from 0 to 1 to use when choosing
            individuals to contact

        Keyword arguments:
        host_origin -- whether to draw from hosts in the origin population
            (as opposed to vectors) (Boolean)
        host_target -- whether to draw from hosts in the target population
            (as opposed to vectors) (Boolean)
        """

        if host_origin:
            index_host,rand = self.getWeightedRandom(
                rand, np.multiply(
                    self.coefficients_hosts[:,self.POPULATION_CONTACT],
                    self.coefficients_hosts[:,self.INFECTED]
                    )
                )
            origin = self.hosts[index_host]
        else:
            index_vector,rand = self.getWeightedRandom(
                rand, np.multiply(
                    self.coefficients_vectors[:,self.POPULATION_CONTACT],
                    self.coefficients_vectors[:,self.INFECTED]
                    )
                )
            origin = self.vectors[index_vector]

        if host_target:
            index_host,rand = target_pop.getWeightedRandom(
                rand,
                target_pop.coefficients_hosts[
                    :,target_pop.RECEIVE_POPULATION_CONTACT
                    ]
                )
            origin.infectHost(target_pop.hosts[index_host])
        else:
            index_vector,rand = target_pop.getWeightedRandom(
                rand,
                target_pop.coefficients_vectors[
                    :,target_pop.RECEIVE_POPULATION_CONTACT
                    ]
                )
            origin.infectVector(target_pop.vectors[index_vector])

    def contactHostHost(self, rand):
        """Contact any two (weighted) random hosts in population.

        Carries out possible infection events from the first organism into the
        second.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individuals to contact

        Returns:
        whether or not the model has changed state (Boolean)
        """

        index_host,rand = self.getWeightedRandom(
            rand, np.multiply(
                self.coefficients_hosts[:,self.CONTACT],
                self.coefficients_hosts[:,self.INFECTED]
                )
            )
        index_other_host,rand = self.getWeightedRandom(
            rand, self.coefficients_hosts[:,self.RECEIVE_CONTACT]
            )

        changed = self.hosts[index_host].infectHost(
            self.hosts[index_other_host]
            )

        return changed

    def contactHostVector(self, rand):
        """Contact a (weighted) random host and vector in population.

        Carries out possible infection events from the first organism into the
        second.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individuals to contact

        Returns:
        whether or not the model has changed state (Boolean)
        """

        index_host,rand = self.getWeightedRandom(
            rand, np.multiply(
                self.coefficients_hosts[:,self.CONTACT],
                self.coefficients_hosts[:,self.INFECTED]
                )
            )
        index_vector,rand = self.getWeightedRandom(
            rand, self.coefficients_vectors[:,self.RECEIVE_CONTACT]
            )
        changed = self.hosts[index_host].infectVector(
            self.vectors[index_vector]
            )

        return changed

    def contactVectorHost(self, rand):
        """Contact a (weighted) random vector and host in population.

        Carries out possible infection events from the first organism into the
        second.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individuals to contact

        Returns:
        whether or not the model has changed state (Boolean)
        """

        index_vector,rand = self.getWeightedRandom(
            rand, np.multiply(
                self.coefficients_vectors[:,self.CONTACT],
                self.coefficients_vectors[:,self.INFECTED]
                )
            )
        index_host,rand = self.getWeightedRandom(
            rand,self.coefficients_hosts[:,self.RECEIVE_CONTACT]
            )

        changed = self.vectors[index_vector].infectHost(
            self.hosts[index_host]
            )

        return changed

    def recoverHost(self, rand):
        """Remove all infections from host at this index.

        If model is protecting upon recovery, add protecion sequence as defined
        by the indexes in the corresponding model parameter. Remove from
        population infected list and add to healthy list.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual to recover
        """

        index_host,rand = self.getWeightedRandom(
            rand,self.coefficients_hosts[:,self.RECOVERY]
            )

        self.hosts[index_host].recover()

    def recoverVector(self, rand):
        """Remove all infections from vector at this index.

        If model is protecting upon recovery, add protecion sequence as defined
        by the indexes in the corresponding model parameter. Remove from
        population infected list and add to healthy list.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual to recover
        """

        index_vector,rand = self.getWeightedRandom(
            rand,self.coefficients_vectors[:,self.RECOVERY]
            )

        self.vectors[index_vector].recover()

    def killHost(self, rand):
        """Add host at this index to dead list, remove it from alive ones.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual to kill
        """

        index_host,rand = self.getWeightedRandom(
            rand, self.lethality_rate_host * self.recovery_rate_host
                * self.coefficients_hosts[:,self.LETHALITY]
            )

        self.dead_hosts.append(self.hosts[index_host])
        self.hosts[index_host].die()

    def killVector(self, rand):
        """Add host at this index to dead list, remove it from alive ones.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual to kill
        """

        index_vector,rand = self.getWeightedRandom(
            rand, self.lethality_rate_vector * self.recovery_rate_vector
                * self.coefficients_vectors[:,self.LETHALITY]
            )

        self.dead_vectors.append(self.vectors[index_vector])
        self.vectors[index_vector].die()

    def dieHost(self, rand):
        """Remove host at this index from alive lists.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual to kill
        """

        index_host,rand = self.getWeightedRandom(
            rand, np.ones( self.coefficients_hosts.shape[0] )
            )

        self.hosts[index_host].die()

    def dieVector(self, rand):
        """Remove vector at this index from alive lists.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual to kill
        """

        index_vector,rand = self.getWeightedRandom(
            rand, np.ones( self.coefficients_vectors.shape[0] )
            )

        self.vectors[index_vector].die()

    def birthHost(self, rand):
        """Add host at this index to population, remove it from alive ones.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual to make parent
        """

        index_host,rand = self.getWeightedRandom(
            rand,self.coefficients_hosts[:,self.NATALITY]
            )

        self.hosts[index_host].birth(rand)

    def birthVector(self, rand):
        """Add host at this index to population, remove it from alive ones.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual to make parent
        """

        index_vector,rand = self.getWeightedRandom(
            rand,self.coefficients_vectors[:,self.NATALITY]
            )

        self.vectors[index_vector].birth(rand)

    def mutateHost(self, rand):
        """Mutate a single, random locus in a random pathogen in the given host.

        Creates a new genotype from a de novo mutation event in the host given.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual in which to choose a pathogen to mutate
        """

        index_host,rand = self.getWeightedRandom(
            rand,self.coefficients_hosts[:,self.MUTATION]
            )
        host = self.hosts[index_host]

        host.mutate(rand)

    def mutateVector(self, rand):
        """Mutate a single, random locus in a random pathogen in given vector.

        Creates a new genotype from a de novo mutation event in the vector
        given.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual in which to choose a pathogen to mutate
        """

        index_vector,rand = self.getWeightedRandom(
            rand,self.coefficients_vectors[:,self.MUTATION]
            )
        vector = self.vectors[index_vector]

        vector.mutate(rand)

    def recombineHost(self, rand):
        """Recombine 2 random pathogen genomes at random locus in given host.

        Creates a new genotype from two random possible pathogens in the host
        given.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual in which to choose pathogens to recombine
        """

        index_host,rand = self.getWeightedRandom(
            rand,self.coefficients_hosts[:,self.MUTATION]
            )
        host = self.hosts[index_host]

        host.recombine(rand)

    def recombineVector(self, rand):
        """Recombine 2 random pathogen genomes at random locus in given vector.

        Creates a new genotype from two random possible pathogens in the vector
        given.

        Arguments:
        rand -- uniform random number from 0 to 1 to use when choosing
            individual in which to choose pathogens to recombine
        """

        index_vector,rand = self.getWeightedRandom(
            rand,self.coefficients_vectors[:,self.MUTATION]
            )
        vector = self.vectors[index_vector]

        vector.recombine(rand)

    def updateHostCoefficients(self):
        """Updates event coefficient values in population's hosts."""
        self.coefficients_hosts = np.zeros( self.coefficients_hosts.shape )
        self.coefficients_hosts[ 1:, self.RECEIVE_CONTACT ] = 1
        self.coefficients_hosts[ 1:, self.RECEIVE_POPULATION_CONTACT ] = 1
        self.coefficients_hosts[ 1:, self.NATALITY ] = 1
        self.coefficients_hosts[ 1:, self.MIGRATION ] = 1

        for h in self.hosts:
            genomes = h.pathogens.keys()
            h.pathogens = {}
            h.sum_fitness = 0
            for g in genomes:
                h.acquirePathogen(g)

    def updateVectorCoefficients(self):
        """Updates event coefficient values in population's vectors."""
        self.coefficients_vectors = np.zeros( self.coefficients_vectors.shape )
        self.coefficients_vectors[ 1:, self.RECEIVE_CONTACT ] = 1
        self.coefficients_vectors[ 1:, self.RECEIVE_POPULATION_CONTACT ] = 1
        self.coefficients_hosts[ 1:, self.NATALITY ] = 1
        self.coefficients_hosts[ 1:, self.MIGRATION ] = 1

        for v in self.vectors:
            genomes = v.pathogens.keys()
            v.pathogens = {}
            v.sum_fitness = 0
            for g in genomes:
                v.acquirePathogen(g)

    def healthyCoefficientRow(self):
        """Returns coefficient values corresponding to a healthy host/vector."""

        v = np.zeros((1,self.NUM_COEFFICIENTS))
        v[ 0, self.RECEIVE_CONTACT ] = 1
        v[ 0, self.RECEIVE_POPULATION_CONTACT ] = 1
        v[ 0, self.NATALITY ] = 1
        v[ 0, self.MIGRATION ] = 1

        return v

    def getWeightedRandom(self, rand, r):
        """Returns index of element chosen from weights and given random number.

        Since sampling from coefficient arrays which contain a dummy first row,
        index is decreased by 1.

        Arguments:
        rand -- 0-1 random number (number)
        r -- array with weights (numpy vector)

        Returns:
        new 0-1 random number (number)
        """

        r_tot = np.sum( r )
        u = rand * r_tot # random uniform number between 0 and total rate
        r_cum = 0
        for i,e in enumerate(r): # for every possible event,
            r_cum += e # add this event's rate to cumulative rate
            if u < r_cum: # if random number is under cumulative rate
                return i-1, ( ( u - r_cum + e ) / e )
