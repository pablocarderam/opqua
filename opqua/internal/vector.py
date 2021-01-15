
"""Contains class Vector."""

import numpy as np
import copy as cp

class Vector(object):
    """Class defines vector entities to be infected by pathogens in model.

    These can infect hosts, the main entities in the model.

    Methods:
    infectHost -- infects given host with a sample of this vector's pathogens
    recover -- removes all infections
    die -- kills this vector
    applyTreatment -- removes all infections with genotypes susceptible to given
        treatment
    """

    def __init__(self, population, id):
        """Create a new Vector.

        Arguments:
        population -- the population this vector belongs to (Population)
        id -- unique identifier for this vector within population (String)
        """
        super(Vector, self).__init__()
        self.id = id
        self.pathogens = {} # Dictionary with all current infections in this
            # vector, with keys=genome strings, values=fitness numbers
        self.population = population
        self.sum_fitness = 0 # sum of all pathogen fitnesses within this vector
        self.protection_sequences = [] # A list of strings this vector is immune
            # to. If a pathogen's genome contains one of these values, it cannot
            # infect this vector.

    def infectHost(self, host):
        """Infect given host with a sample of this vector's pathogens.

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
        """Remove all infections from this vector.

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
        self.sum_fitness = 0
        if self in self.population.infected_vectors:
            self.population.infected_vectors.remove(self)

    def die(self):
        """Add vector to population's dead list, remove it from alive ones."""

        self.population.dead_vectors.append(self)
        self.population.vectors.remove(self)

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
            self.sum_fitness -= self.pathogens[g]
            del self.pathogens[g]

        if ( len(self.pathogens) == 0
                and self in self.population.infected_vectors ):
            self.population.infected_vectors.remove(self)
