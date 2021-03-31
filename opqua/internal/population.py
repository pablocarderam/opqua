
"""Contains class Population."""

import numpy as np
import copy as cp
import random

from opqua.internal.host import Host
from opqua.internal.vector import Vector

class Population(object):
    """Class defines a population with hosts, vectors, and specific parameters.

    Methods:
    setSetup -- assigns a given set of parameters to this population
    copyState -- returns a slimmed-down version of the current population state
    addHosts -- adds hosts to the population
    addVectors -- adds vectors to the population
    newHostGroup -- returns a list of random (healthy or any) hosts
    newVectorGroup -- returns a list of random (healthy or any) vectors
    removeHosts -- removes hosts from the population
    removeVectors -- removes vectors from the population
    addPathogensToHosts -- adds pathogens with specified genomes to hosts
    addPathogensToVectors -- adds pathogens with specified genomes to vectors
    recoverHost --removes all infections from given host
    recoverVector -- removes all infections from given vector
    treatHosts -- removes infections susceptible to given treatment from hosts
    treatVectors -- removes infections susceptible to treatment from vectors
    protectHosts -- adds protection sequence to hosts
    protectVectors -- adds protection sequence to vectors
    wipeProtectionHosts -- removes all protection sequences from hosts
    wipeProtectionVectors -- removes all protection sequences from hosts
    setNeighbor -- sets migration rate from this population towards another
    migrate -- transfers hosts and/or vectors from this population to a neighbor
    contactInfectedHostAnyHost -- carries out a contact event between a random
        infected host and any random host in population
    contactInfectedHostAnyVector -- carries out a contact event between a random
        infected host and any random vector in population
    contactHealthyHostInfectedVector -- carries out a contact event between a
        random healthy host and a random infected vector in population
    mutateHost -- mutates a single locus in a random pathogen in a host
    mutateVector -- mutates a single locus in a random pathogen in a vector
    recombineHost -- recombines two random pathogens in a host
    recombineVector -- recombines two random pathogens in a host
    updateHostFitness -- updates fitness of pathogens in population's hosts
    updateVectorFitness -- updates fitness of pathogens in population's vectors
    """

    def __init__(self, id, setup, num_hosts, num_vectors, slim=False):
        """Create a new Population.

        Arguments:
        id -- unique identifier for this population in the model (String)
        setup -- setup object with parameters for this population (Setup)
        num_hosts -- number of hosts to initialize population with (int)
        num_vectors -- number of hosts to initialize population with (int)
        slim -- whether to create a slimmed-down representation of the
            population for data storage (only ID, host and vector lists)
            (Boolean, default False)
        """
        super(Population, self).__init__()

        self.id = id

        if not slim:
                # if not slimmed down for data storage, save other attributes
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
            self.infected_hosts = []
                # contains live, infected hosts (also in hosts list)
            self.healthy_hosts = self.hosts.copy()
                # contains live, healthy hosts (also in hosts list)
            self.infected_vectors = []
                # contains live, infected vectors (also in vectors list)
            self.dead_hosts = [] # contains dead hosts (not in hosts list)
            self.dead_vectors = [] # contains dead vectors (not in vectors list)
            self.neighbors = {} # dictionary with neighboring populations,
                # keys=Population, values=migration rate from this population to
                # neighboring population

            self.total_migration_rate = 0 # sum of all migration rates from this
                # population to neighbors

            self.setSetup(setup)

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

        copy = Population(self.id, None, 0, 0, slim=True)
        if host_sampling > 0:
            host_sample = random.sample(
                self.hosts, int( len(self.hosts) / host_sampling )
                )
            dead_host_sample = random.sample(
                self.dead_hosts, int( len(self.dead_hosts) / host_sampling )
                )
            copy.hosts = [ h.copyState() for h in host_sample ]
            copy.dead_hosts = [ h.copyState() for h in dead_host_sample ]
        else:
            copy.hosts = [ h.copyState() for h in self.hosts ]
            copy.dead_hosts = [ h.copyState() for h in self.dead_hosts ]

        if vector_sampling > 0:
            vector_sample = random.sample(
                self.vectors, int( len(self.vectors) / vector_sampling )
                )
            dead_vector_sample = random.sample(
                self.dead_vectors, int( len(self.dead_vectors)/vector_sampling )
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
        self.updateHostFitness()
        self.fitnessVector = setup.fitnessVector
        self.updateVectorFitness()
        self.lethalityHost = setup.lethalityHost
        self.updateHostLethality()
        self.lethalityVector = setup.lethalityVector
        self.updateVectorLethality()

        self.contact_rate_host_vector = setup.contact_rate_host_vector
        self.contact_rate_host_host = setup.contact_rate_host_host
            # contact rate assumes fixed area--large populations are dense
            # populations, so contact scales linearly with both host and vector
            # populations. If you don't want this to happen, modify the
            # population's contact rate accordingly.
        self.mean_inoculum_host = setup.mean_inoculum_host
        self.mean_inoculum_vector = setup.mean_inoculum_vector
        self.recovery_rate_host = setup.recovery_rate_host
        self.recovery_rate_vector = setup.recovery_rate_vector
        self.recombine_in_host = setup.recombine_in_host
        self.recombine_in_vector = setup.recombine_in_vector
        self.mutate_in_host = setup.mutate_in_host
        self.mutate_in_vector = setup.mutate_in_vector
        self.death_rate_host = setup.death_rate_host
        self.death_rate_vector = setup.death_rate_vector
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

    def newHostGroup(self, num_hosts, healthy=True):
        """Return a list of random (healthy or any) hosts in population.

        Arguments:
        num_vectors -- number of vectors to be sampled randomly (int)

        Keyword arguments:
        healthy -- whether to sample healthy hosts only (default True; Boolean)

        Returns:
        list containing sampled hosts
        """

        possible_hosts = []
        if healthy:
            if len(self.healthy_hosts) >= num_hosts:
                possible_hosts = self.healthy_hosts
            else:
                raise ValueError(
                    "You're asking for " + str(num_hosts)
                    + " healthy hosts, but population " + str(self.id)
                    + " only has " + str( len(self.healthy_hosts) ) + "."
                    )
        else:
            possible_hosts = self.hosts

        hosts = np.random.choice(possible_hosts, num_hosts, replace=False)

        return hosts

    def newVectorGroup(self, num_vectors, healthy=True):
        """Return a list of random (healthy or any) vectors in population.

        Arguments:
        num_vectors -- number of vectors to be sampled randomly (int)

        Keyword arguments:
        healthy -- whether to sample healthy vectors only (default True;
            Boolean)

        Returns:
        list containing sampled vectors
        """

        possible_vectors = []
        if healthy:
            if len(self.vectors) - len(self.infected_vectors) >= num_vectors:
                for vector in self.vectors:
                    if vector not in self.infected_vectors:
                        possible_vectors.append(vector)
            else:
                raise ValueError(
                    "You're asking for " + str(num_vectors)
                    + " healthy vectors, but population " + str(self.id)
                    + " only has "
                    + str( len(self.vectors) - len(self.infected_vectors) )
                    + "."
                    )
        else:
            possible_vectors = self.vectors

        vectors = np.random.choice(possible_vectors, num_vectors, replace=False)

        return vectors

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
        else:
            for _ in range(num_hosts_or_list):
                host_removed = np.random.choice(self.hosts)
                if host_removed in self.infected_hosts:
                    self.infected_hosts.remove( host_removed )
                else:
                    self.healthy_hosts.remove( host_removed )

                self.hosts.remove( host_removed )

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
        else:
            for _ in range(num_vectors):
                vector_removed = np.random.choice(self.vectors)
                if vector_removed in self.infected_vectors:
                    self.infected_vectors.remove( vector_removed )

                self.vectors.remove( vector_removed )

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

        for genome in genomes_numbers:
            if len(genome) == self.num_loci and all( [
                allele in self.possible_alleles[i]
                    for i,allele in enumerate(genome)
                ] ):
                new_fitness = self.fitnessHost(genome)
                if len(hosts) == 0:
                    hosts = self.hosts

                rand_hosts = np.random.choice(
                    hosts, int(genomes_numbers[genome]), replace=False
                    )
                for rand_host in rand_hosts:
                    if rand_host not in self.infected_hosts:
                        self.infected_hosts.append(rand_host)
                        self.healthy_hosts.remove(rand_host)

                    if genome not in rand_host.pathogens:
                        rand_host.pathogens[genome] = new_fitness
                        rand_host.sum_fitness += new_fitness
                        if rand_host.pathogens[genome] > rand_host.max_fitness:
                            rand_host.max_fitness = rand_host.pathogens[genome]

            else:
                raise ValueError('Genome ' + genome + ' must be of length '
                    + str(self.num_loci)
                    + ' and contain only the following characters at each '
                    + 'position: ' + self.possible_alleles + ' .')

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

        for genome in genomes_numbers:
            if len(genome) == self.num_loci and all( [
                allele in self.possible_alleles[i]
                    for i,allele in enumerate(genome)
                ] ):
                new_fitness = self.fitnessVector(genome)
                if len(vectors) == 0:
                    vectors = self.vectors

                rand_vectors = np.random.choice(
                    vectors, int(genomes_numbers[genome]), replace=False
                    )
                for rand_vector in rand_vectors:
                    if rand_vector not in self.infected_vectors:
                        self.infected_vectors.append(rand_vector)

                    if genome not in rand_vector.pathogens:
                        rand_vector.pathogens[genome] = new_fitness
                        rand_vector.sum_fitness += new_fitness
                        if rand_vector.pathogens[genome] > rand_vector.max_fitness:
                            rand_vector.max_fitness = rand_vector.pathogens[genome]
            else:
                raise ValueError('Genome ' + genome + ' must be of length '
                    + str(self.num_loci)
                    + ' and contain only the following characters at each '
                    + 'position: ' + self.possible_alleles + ' .')

    def recoverHost(self, index_host):
        """Remove all infections from host at this index.

        If model is protecting upon recovery, add protecion sequence as defined
        by the indexes in the corresponding model parameter. Remove from
        population infected list and add to healthy list.

        Arguments:
        index_host -- index of host in infected_hosts (int)
        """

        self.infected_hosts[index_host].recover()

    def recoverVector(self, index_vector):
        """Remove all infections from vector at this index.

        If model is protecting upon recovery, add protecion sequence as defined
        by the indexes in the corresponding model parameter. Remove from
        population infected list and add to healthy list.

        Arguments:
        index_vector -- index of vector in infected_vectors (int)
        """

        self.infected_vectors[index_vector].recover()

    def killHost(self, index_host, rand):
        """Add host at this index to dead list, remove it from alive ones.

        Arguments:
        index_host -- index of host in infected_hosts (int)
        rand -- a random number between 0 and 1 used to determine death
        """

        if infected_hosts[index_host].max_lethality > rand:
            self.infected_hosts[index_host].die()

    def killVector(self, index_vector, rand):
        """Add host at this index to dead list, remove it from alive ones.

        Arguments:
        index_vector -- index of vector in infected_vectors (int)
        rand -- a random number between 0 and 1 used to determine death
        """

        if infected_hosts[index_host].max_lethality > rand:
            self.infected_vectors[index_vector].die()

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
            possible_infected_vectors[vector].applyTreatment(resistance_seqs)

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

    def setNeighbor(self, neighbor, rate):
        """Set migration rate from this population towards another one.

         Arguments:
         neighbor -- population towards which migration rate will be specified
            (Population)
         rate -- migration rate from this population to the neighbor (number)
         """

        if neighbor in self.neighbors:
            self.total_migration_rate -= self.neighbors[neighbor]

        self.neighbors[neighbor] = rate
        self.total_migration_rate += rate

    def migrate(self, target_pop, num_hosts, num_vectors):
        """Transfer hosts and/or vectors from this population to another.

        Arguments:
        target_pop -- population towards which migration will occur (Population)
        num_hosts -- number of hosts to transfer (int)
        num_vectors -- number of vectors to transfer (int)
        """

        migrating_hosts = np.random.choice(self.hosts, num_hosts, replace=False)
        for host in migrating_hosts:
            self.hosts.remove(host)
            target_pop.hosts.append(host)
            host.population = target_pop
            if host in self.infected_hosts:
                self.infected_hosts.remove(host)
                target_pop.infected_hosts.append(host)
            else:
                self.healthy_hosts.remove(host)
                target_pop.healthy_hosts.append(host)

        migrating_vectors = np.random.choice(
            self.vectors,num_vectors, replace=False
            )
        for vector in migrating_vectors:
            self.vectors.remove(vector)
            target_pop.vectors.append(vector)
            vector.population = target_pop
            if vector in self.infected_vectors:
                self.infected_vectors.remove(vector)
                target_pop.infected_vectors.append(vector)

    def contactInfectedHostAnyHost(self, index_infected_host, index_other_host):
        """Contact a random infected host and any random host in population.

        Carries out possible infection events from each organism into the other.

        Arguments:
        index_infected_host -- index of infected host in infected_hosts (int)
        index_other_host -- index of second host in hosts (int)

        Returns:
        whether or not the model has changed state (Boolean)
        """
        changed = self.infected_hosts[index_infected_host].infectHost(self.hosts[index_other_host])

        return changed

    def contactInfectedHostAnyVector(self, index_host, index_vector):
        """Contact a random infected host and any random vector in population.

        Carries out possible infection events from each organism into the other.

        Arguments:
        index_host -- index of infected host in infected_hosts (int)
        index_vector -- index of vector in vectors (int)

        Returns:
        whether or not the model has changed state (Boolean)
        """

        changed = self.infected_hosts[index_host].infectVector(self.vectors[index_vector])

        return changed

    def contactHealthyHostInfectedVector(self, index_host, index_vector):
        """Contact a random healthy host and a random infected vector.

        Carries out possible infection events from vector into host.

        Arguments:
        index_host -- index of host in healthy_hosts (int)
        index_vector -- index of infected vector in infected_vectors (int)

        Returns:
        whether or not the model has changed state (Boolean)
        """

        changed = self.infected_vectors[index_vector].infectHost(
            self.healthy_hosts[index_host]
            )

        return changed

    def mutateHost(self, index_host):
        """Mutate a single, random locus in a random pathogen in the given host.

        Creates a new genotype from a de novo mutation event in the host given.

        Arguments:
        index_host -- index of host in infected_hosts (int)
        """

        host = self.infected_hosts[index_host]

        if len(host.pathogens) > 0:
            old_genome = np.random.choice( list( host.pathogens.keys() ) )
            mut_index = np.random.randint( self.num_loci )
            new_genome = old_genome[0:mut_index] + np.random.choice(
                list(self.possible_alleles[mut_index])
                ) + old_genome[mut_index+1:]
            if new_genome not in host.pathogens:
                host.pathogens[new_genome] = self.fitnessHost(new_genome)
                host.sum_fitness += host.pathogens[new_genome]
                if host.pathogens[new_genome] > host.max_fitness:
                    host.max_fitness = host.pathogens[new_genome]

                lethality = self.lethalityHost(new_genome)
                if lethality > host.max_lethality:
                    host.max_lethality = lethality

    def mutateVector(self, index_vector):
        """Mutate a single, random locus in a random pathogen in given vector.

        Creates a new genotype from a de novo mutation event in the vector
        given.

        Arguments:
        index_vector -- index of vector in infected_vectors (int)
        """

        vector = self.infected_vectors[index_vector]

        if len(vector.pathogens) > 0:
            old_genome = np.random.choice( list( vector.pathogens.keys() ) )
            mut_index = np.random.randint( self.num_loci )
            new_genome = old_genome[0:mut_index] + np.random.choice(
                list(self.possible_alleles[mut_index])
                ) + old_genome[mut_index+1:]
            if new_genome not in vector.pathogens:
                vector.pathogens[new_genome] = self.fitnessHost(new_genome)
                vector.sum_fitness += vector.pathogens[new_genome]
                if vector.pathogens[new_genome] > vector.max_fitness:
                    vector.max_fitness = vector.pathogens[new_genome]

                lethality = self.lethalityVector(new_genome)
                if lethality > vector.max_lethality:
                    vector.max_lethality = lethality

    def recombineHost(self, index_host):
        """Recombine two random pathogen genomes at random locus in given host.

        Creates a new genotype from two random possible pathogens in the host
        given.

        Arguments:
        index_host -- index of host in infected_hosts (int)
        """

        host = self.infected_hosts[index_host]

        if len(host.pathogens) > 1:
            parents = np.random.choice( list( host.pathogens.keys() ), 2,
                replace=False )
            recom_index = np.random.randint( self.num_loci )
            new_genome = parents[0][0:recom_index] + parents[1][recom_index:]
            if new_genome not in host.pathogens:
                host.pathogens[new_genome] = self.fitnessHost(new_genome)
                host.sum_fitness += host.pathogens[new_genome]
                if host.pathogens[new_genome] > host.max_fitness:
                    host.max_fitness = host.pathogens[new_genome]

                lethality = self.lethalityHost(new_genome)
                if lethality > host.max_lethality:
                    host.max_lethality = lethality

    def recombineVector(self, index_vector):
        """Recombine 2 random pathogen genomes at random locus in given vector.

        Creates a new genotype from two random possible pathogens in the vector
        given.

        Arguments:
        index_vector -- index of vector in infected_vectors (int)
        """

        vector = self.infected_vectors[index_vector]

        if len(vector.pathogens) > 1:
            parents = np.random.choice( list( vector.pathogens.keys() ), 2,
                replace=False )
            recom_index = np.random.randint( self.num_loci )
            new_genome = parents[0][0:recom_index] + parents[1][recom_index:]
            if new_genome not in vector.pathogens:
                vector.pathogens[new_genome] = self.fitnessHost(new_genome)
                vector.sum_fitness += vector.pathogens[new_genome]
                if vector.pathogens[new_genome] > vector.max_fitness:
                    vector.max_fitness = vector.pathogens[new_genome]

                lethality = self.lethalityVector(new_genome)
                if lethality > vector.max_lethality:
                    vector.max_lethality = lethality

    def updateHostFitness(self):
        """Updates fitness values of all pathogens in hosts for this population.

        """
        for h in self.hosts:
            h.pathogens = {
                g:self.fitnessHost(g) for g,f in h.pathogens.items()
                }
            fitnesses = np.array( list( h.pathogens.values() ) )
            if len(fitnesses) > 0:
                h.sum_fitness = fitnesses.sum()
                h.max_fitness = fitnesses.max()
            else:
                h.sum_fitness = 0
                h.max_fitness = 0

    def updateVectorFitness(self):
        """Updates fitness values of all pathogens in vectors for this
        population.

        """
        for v in self.vectors:
            v.pathogens = {
                g:self.fitnessVector(g) for g,f in v.pathogens.items()
                }
            fitnesses = np.array( list( v.pathogens.values() ) )
            if len(fitnesses) > 0:
                v.sum_fitness = fitnesses.sum()
                v.max_fitness = fitnesses.max()
            else:
                v.sum_fitness = 0
                v.max_fitness = 0

    def updateHostLethality(self):
        """Updates lethality values of all pathogens in hosts for population.

        """
        for h in self.hosts:
            if len(h.pathogens) > 0:
                h.max_lethality = max(
                    [ self.lethalityHost(g) for g,f in h.pathogens.items() ]
                    )
            else:
                h.max_lethality = 0

    def updateVectorLethality(self):
        """Updates lethality values of all pathogens in vectors for this
        population.

        """
        for v in self.vectors:
            if len(v.pathogens) > 0:
                v.max_lethality = max(
                    [ self.lethalityVector(g) for g,f in v.pathogens.items() ]
                    )
            else:
                v.max_lethality = 0
