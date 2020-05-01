# TODO: Pathogen genome influences transmission probability
"""Contains class Host."""

import numpy as np
import copy as cp

class Host(object):
    """Class defines main entities to be infected by pathogens in model.

    Methods:
    infectVector -- infects given vector with a sample of this host's pathogens
    infectHost -- infects given host with a sample of this host's pathogens
    recover -- removes all infections
    die -- kills this host
    applyTreatment -- removes all infections with genotypes susceptible to given
        treatment
    """

    def __init__(self, population, id):
        """Create a new Host.

        Arguments:
        population -- the population this host belongs to (Population)
        id -- unique identifier for this host within population (String)
        """
        super(Host, self).__init__()
        self.id = id
        self.pathogens = {} # Dictionary with all current infections in this
            # host, with keys=genome strings, values=fitness numbers
        self.population = population
        self.sum_fitness = 0 # sum of all pathogen fitnesses within this host
        self.protection_sequences = [] # A list of strings this host is immune
            # to. If a pathogen's genome contains one of these values, it cannot
            # infect this host.

    def infectVector(self, vector):
        """Infect given vector with a sample of this host's pathogens.

        Each pathogen in the infector is sampled as present or absent in the
        inoculum by drawing from a Poisson distribution with a mean equal to the
        mean inoculum size of the organism being infected weighed by each
        genome's fitness as a fraction of the total in the infector) as the
        probability of each trial. Each pathogen present in the inoculum will be
        added to the infected organism, if it does not have protection from the
        pathogen's genome. Fitnesses are computed for the pathogens' genomes in
        the infected organism, and the organism is included in the poplation's
        infected list if appropriate.

        Arguments:
        vector -- the vector to be infected (Vector)

        Returns:
        whether or not the model has changed state (Boolean)
        """

        changed = False
        for genome,fitness in self.pathogens.items():
            if ( genome not in vector.pathogens.keys()
                    and not any(
                        [ p in genome for p in vector.protection_sequences ]
                    )
                    and np.random.poisson(
                        self.population.mean_inoculum_vector
                            * fitness / self.sum_fitness, 1
                        )
                    > 0 ):
                vector.pathogens[genome] = vector.population.fitnessVector(
                    genome
                    )
                vector.sum_fitness += vector.pathogens[genome]
                changed = True
                if vector not in vector.population.infected_vectors:
                    vector.population.infected_vectors.append(vector)

        return changed

    def infectHost(self, host):
        """Infect given host with a sample of this host's pathogens.

        Each pathogen in the infector is sampled as present or absent in the
        inoculum by drawing from a Poisson distribution with a mean equal to the
        mean inoculum size of the organism being infected weighed by each
        genome's fitness as a fraction of the total in the infector) as the
        probability of each trial. Each pathogen present in the inoculum will be
        added to the infected organism, if it does not have protection from the
        pathogen's genome. Fitnesses are computed for the pathogens' genomes in
        the infected organism, and the organism is included in the poplation's
        infected list if appropriate.

        Arguments:
        host -- the host to be infected (Host)

        Returns:
        whether or not the model has changed state (Boolean)
        """

        changed = False
        for genome,fitness in self.pathogens.items():
            if ( genome not in host.pathogens.keys()
                    and not any(
                        [ p in genome for p in host.protection_sequences ]
                    )
                    and np.random.poisson(
                        self.population.mean_inoculum_host
                            * fitness / self.sum_fitness, 1
                        )
                    > 0 ):
                host.pathogens[genome] = host.population.fitnessHost(genome)
                host.sum_fitness += host.pathogens[genome]
                changed = True
                if host not in host.population.infected_hosts:
                    host.population.infected_hosts.append(host)
                    host.population.healthy_hosts.remove(host)

        return changed

    def recover(self):
        """Remove all infections from this host.

        If model is protecting upon recovery, add protecion sequence as defined
        by the indexes in the corresponding model parameter. Remove from
        population infected list and add to healthy list.
        """

        if self.population.protection_upon_recovery_host:
            for genome in self.pathogens:
                seq = genome[self.protection_upon_recovery_host[0]
                    :self.protection_upon_recovery_host[1]]
                if seq not in self.protection_sequences:
                    self.protection_sequences.append(seq)

        self.pathogens = {}
        if self in self.population.infected_hosts:
            self.population.infected_hosts.remove(self)
            self.population.healthy_hosts.append(self)

    def die(self):
        """Add host to population's dead list, remove it from alive ones."""

        self.population.dead_hosts.append(self)
        if self in self.population.infected_hosts:
            self.population.infected_hosts.remove(host)
        else:
            self.population.healthy_hosts.remove(self)

        self.population.hosts.remove(self)

    def applyTreatment(self, resistance_seqs):
        """Remove all infections with genotypes susceptible to given treatment.

        Pathogens are removed if they are missing at least one of the sequences
        in resistance_seqs from their genome. Removes this organism from
        population infected list and adds to healthy list if appropriate.

        Arguments:
        resistance_seqs -- contains sequences required for treatment resistance
            (list of Strings)
        """

        genomes_removed = []
        for genome in self.pathogens:
            for seq in resistance_seqs:
                if seq not in genome:
                    genomes_removed.append(genome)
                    break

        for g in genomes_removed:
            del self.pathogens[g]

        if len(self.pathogens) == 0 and self in self.population.infected_hosts:
            self.population.infected_hosts.remove(self)
            self.population.healthy_hosts.append(self)
