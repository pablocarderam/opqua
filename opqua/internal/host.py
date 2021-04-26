
"""Contains class Host."""

import numpy as np
import copy as cp

class Host(object):
    """Class defines main entities to be infected by pathogens in model.

    Methods:
    copyState -- returns a slimmed-down version of the current host state
    infectVector -- infects given vector with a sample of this host's pathogens
    infectHost -- infects given host with a sample of this host's pathogens
    recover -- removes all infections
    die -- kills this host
    applyTreatment -- removes all infections with genotypes susceptible to given
        treatment
    """

    def __init__(self, population, id, slim=False):
        """Create a new Host.

        Arguments:
        population -- the population this host belongs to (Population)
        id -- unique identifier for this host within population (String)
        slim -- whether to create a slimmed-down representation of the
            population for data storage (only ID, host and vector lists)
            (Boolean, default False)
        """
        super(Host, self).__init__()
        self.id = id

        if not slim:
                # if not slimmed down for data storage, save other attributes
            self.pathogens = {} # Dictionary with all current infections in this
                # host, with keys=genome strings, values=fitness numbers
            self.protection_sequences = [] # A list of strings this host is
                # immune to. If a pathogen's genome contains one of these
                # values, it cannot infect this host.
            self.population = population
            self.sum_fitness = 0
                # sum of all pathogen fitnesses within this host
            self.coefficient_index = population.coefficients_hosts.shape[0] \
                if population.coefficients_hosts.ndim > 1 else 1
                # index in population's coefficient array, not same as id

            population.coefficients_hosts = np.vstack( (
                population.coefficients_hosts,
                np.zeros( population.NUM_COEFFICIENTS )
                ) ) # adds a row to coefficient array

    def copyState(self):
        """Returns a slimmed-down representation of the current host state.

        Returns:
        Host object with current pathogens and protection_sequences.
        """

        copy = Host(None, self.id, slim=True)
        copy.pathogens = self.pathogens.copy()
        copy.protection_sequences = self.protection_sequences.copy()

        return copy

    def acquirePathogen(self, genome):
        """Adds given genome to this host's pathogens.

        Modifies event coefficient matrix accordingly.

        Arguments:
        genome -- the genome to be added (String)
        """
        self.pathogens[genome] = self.population.fitnessHost(genome)
        old_sum_fitness = self.sum_fitness
        self.sum_fitness += self.pathogens[genome]
        self.population.coefficients_hosts[self.coefficient_index,:] = (
            self.population.coefficients_hosts[
                self.coefficient_index,:
                ]
            * old_sum_fitness / self.sum_fitness ) + ( np.array([
                    # positions dependent on class constants
                0,
                1-self.population.contactHost(genome),
                self.population.lethalityHost(genome),
                self.population.recoveryHost(genome),
                self.population.migrationHost(genome),
                1-self.population.populationContactHost(genome),
                self.population.mutationHost(genome),
                self.population.recombinationHost(genome)
            ]) * self.pathogens[genome] / self.sum_fitness )

        self.population.coefficients_hosts[
            self.coefficient_index,self.population.INFECTED
            ] = 1

        if self not in self.population.infected_hosts:
            self.population.infected_hosts.append(self)
            self.population.healthy_hosts.remove(self)

    def infectHost(self, host):
        """Infect given host with a sample of this host's pathogens.

        Each pathogen in the infector is sampled as present or absent in the
        inoculum by drawing from a Poisson distribution with a mean equal to the
        mean inoculum size of the organism being infected weighted by each
        genome's fitness as a fraction of the total in the infector as the
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
            if genome not in host.pathogens.keys() and not any(
                    [ p in genome for p in host.protection_sequences ]
                    ) and ( np.random.poisson(
                    self.population.mean_inoculum_host
                    * fitness / self.sum_fitness, 1
                    ) > 0 ):
                host.acquirePathogen(genome)
                changed = True

        return changed

    def infectVector(self, vector):
        """Infect given host with a sample of this vector's pathogens.

        Each pathogen in the infector is sampled as present or absent in the
        inoculum by drawing from a Poisson distribution with a mean equal to the
        mean inoculum size of the organism being infected weighted by each
        genome's fitness as a fraction of the total in the infector as the
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
            if genome not in vector.pathogens.keys() and not any(
                    [ p in genome for p in vector.protection_sequences ]
                    ) and ( np.random.poisson(
                    self.population.mean_inoculum_vector
                    * fitness / self.sum_fitness, 1
                    ) > 0 ):
                vector.acquirePathogen(genome)
                changed = True

        return changed

    def recover(self):
        """Remove all infections from this host.

        If model is protecting upon recovery, add protecion sequence as defined
        by the indexes in the corresponding model parameter. Remove from
        population infected list and add to healthy list.
        """

        if self in self.population.infected_hosts:
            if self.population.protection_upon_recovery_host:
                for genome in self.pathogens:
                    seq = genome[self.protection_upon_recovery_host[0]
                        :self.protection_upon_recovery_host[1]]
                    if seq not in self.protection_sequences:
                        self.protection_sequences.append(seq)

            self.pathogens = {}
            self.sum_fitness = 0
            self.population.coefficients_hosts[
                self.coefficient_index,:
                ] = np.zeros( self.population.NUM_COEFFICIENTS )

            self.population.infected_hosts.remove(self)
            self.population.healthy_hosts.append(self)

    def die(self):
        """Add host to population's dead list, remove it from alive ones."""

        self.population.dead_hosts.append(self)
        if self in self.population.infected_hosts:
            self.population.infected_hosts.remove(host)
            self.population.coefficients_hosts[
                self.coefficient_index,:
                ] = np.zeros( self.population.NUM_COEFFICIENTS )
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

        genomes_remaining = []
        for genome in self.pathogens:
            for seq in resistance_seqs:
                if seq in genome:
                    genomes_remaining.append(genome)
                    break

        if len(genomes_remaining) == 0:
            self.recover()
        else:
            self.pathogens = {}
            self.sum_fitness = 0
            self.population.coefficients_hosts[
                self.coefficient_index,:
                ] = np.zeros( self.population.NUM_COEFFICIENTS )
            for genome in genomes_remaining:
                self.acquirePathogen(genome)

    def mutate(self,rand):
        """Mutate a single, random locus in a random pathogen.

        Creates a new genotype from a de novo mutation event.
        """

        genomes = list( self.pathogens.keys() )
        weights = [
            self.population.mutationHost(g)
            * self.population.fitnessHost(g) for g in genomes
            ]
        index_genome,rand = self.population.getWeightedRandom(
            rand,weights
            )

        old_genome = genomes[index_genome]
        mut_index = int( rand * self.population.num_loci )
        if old_genome[mut_index] != '/':
            new_genome = old_genome[0:mut_index] + np.random.choice(
                list(self.population.possible_alleles[mut_index])
                ) + old_genome[mut_index+1:]
            if new_genome not in self.pathogens:
                self.acquirePathogen(new_genome)

    def recombine(self,rand):
        """Recombine two random pathogen genomes at random locus.

        Creates a new genotype from two random possible pathogens.
        """

        genomes = list( self.pathogens.keys() )
        weights = [
            self.population.recombinationHost(g)
            * self.population.fitnessHost(g) for g in genomes
            ]
        index_genome,rand = self.population.getWeightedRandom(
            rand,weights
            )
        index_other_genome,rand = self.population.getWeightedRandom(
            rand,weights
            )

        num_evts = np.random.poisson( self.population.num_crossover_host )
        loci = np.random.randint( 0, self.population.num_loci, num_evts )

        children = [ genomes[index_genome], genomes[index_other_genome] ]

        for l in loci:
            children[0] = children[0][0:l] + children[1][l:]
            children[1] = children[1][0:l] + children[0][l:]

        children = [ genome.split('/') for genome in genomes ]
        parent = np.random.randint( 0, 2, len( children[0] ) )

        children = [
            [
                children[ parent[i] ][i]
                for i in range( len( children[0] ) )
                ].join('/'),
            [
                children[ not parent[i] ][i]
                for i in range( len( children[1] ) )
                ].join('/')
            ]

        for new_genome in children:
            if new_genome not in self.pathogens:
                self.acquirePathogen(new_genome)
