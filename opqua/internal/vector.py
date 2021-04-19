
"""Contains class Vector."""

import numpy as np
import copy as cp

class Vector(object):
    """Class defines vector entities to be infected by pathogens in model.

    These can infect hosts, the main entities in the model.

    Methods:
    copyState -- returns a slimmed-down version of the current vector state
    infectHost -- infects given host with a sample of this vector's pathogens
    recover -- removes all infections
    die -- kills this vector
    applyTreatment -- removes all infections with genotypes susceptible to given
        treatment
    """

    def __init__(self, population, id, slim=False):
        """Create a new Vector.

        Arguments:
        population -- the population this vector belongs to (Population)
        id -- unique identifier for this vector within population (String)
        slim -- whether to create a slimmed-down representation of the
            population for data storage (only ID, host and vector lists)
            (Boolean, default False)
        """
        super(Vector, self).__init__()
        self.id = id

        if not slim:
                # if not slimmed down for data storage, save other attributes
            self.pathogens = {} # Dictionary with all current infections in this
                # vector, with keys=genome strings, values=fitness numbers
            self.protection_sequences = [] # A list of strings this vector is
                # immune to. If a pathogen's genome contains one of these
                # values, it cannot infect this vector.
            self.population = population
            self.sum_fitness = 0
                # sum of all pathogen fitnesses within this vector
            self.coefficient_index = population.coefficients_vectors.shape[0] \
                if population.coefficients_vectors.ndim > 1 else 1
                # index in population's coefficient array, not same as id

            population.coefficients_vectors = np.vstack( (
                population.coefficients_vectors,
                np.zeros( population.NUM_COEFFICIENTS )
                ) ) # adds a row to coefficient array

    def copyState(self):
        """Returns a slimmed-down representation of the current vector state.

        Returns:
        Vector object with current pathogens and protection_sequences.
        """

        copy = Vector(None, self.id, slim=True)
        copy.pathogens = self.pathogens.copy()
        copy.protection_sequences = self.protection_sequences.copy()

        return copy

    def acquirePathogen(self, genome):
        """Adds given genome to this vector's pathogens.

        Modifies event coefficient matrix accordingly.

        Arguments:
        genome -- the genome to be added (String)

        Returns:
        whether or not the model has changed state (Boolean)
        """

        self.pathogens[genome] = self.population.fitnessVector(genome)
        old_sum_fitness = self.sum_fitness
        self.sum_fitness += self.pathogens[genome]
        self.population.coefficients_vectors[self.coefficient_index,:] = (
            self.population.coefficients_vectors[
                self.coefficient_index,:
                ]
            * old_sum_fitness / self.sum_fitness ) + ( np.array([
                    # positions dependent on class constants
                0,
                1-self.population.contactVector(genome),
                self.population.lethalityVector(genome),
                self.population.recoveryVector(genome),
                self.population.migrationVector(genome),
                1-self.population.populationContactVector(genome),
                self.population.mutationVector(genome),
                self.population.recombinationVector(genome)
            ]) * self.pathogens[genome] / self.sum_fitness )

        self.population.coefficients_vectors[
            self.coefficient_index,self.population.INFECTED
            ] = 1

        if self not in self.population.infected_vectors:
            self.population.infected_vectors.append(self)

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
        vector -- the vector to be infected (Vector)

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

    def recover(self):
        """Remove all infections from this vector.

        If model is protecting upon recovery, add protection sequence as defined
        by the indexes in the corresponding model parameter. Remove from
        population infected list and add to healthy list.
        """

        if self in self.population.infected_vectors:
            if self.population.protection_upon_recovery_vector:
                for genome in self.pathogens:
                    seq = genome[self.protection_upon_recovery_vector[0]
                        :self.protection_upon_recovery_vector[1]]
                    if seq not in self.protection_sequences:
                        self.protection_sequences.append(seq)

            self.pathogens = {}
            self.sum_fitness = 0
            self.population.coefficients_vectors[
                self.coefficient_index,:
                ] = np.zeros( self.population.NUM_COEFFICIENTS )

            self.population.infected_vectors.remove(self)

    def die(self):
        """Add vector to population's dead list, remove it from alive ones."""

        self.population.dead_vectors.append(self)
        if self in self.population.infected_vectors:
            self.population.infected_vectors.remove(self)
            self.population.coefficients_vectors[
                self.coefficient_index,:
                ] = np.zeros( self.population.NUM_COEFFICIENTS )

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
            self.population.coefficients_vectors[
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
            self.population.mutationVector(g)
            * self.population.fitnessVector(g) for g in genomes
            ]
        index_genome,rand = self.population.getWeightedRandom(
            rand,weights
            )

        old_genome = genomes[index_genome]
        mut_index = int( rand * self.population.num_loci )
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
            self.population.recombinationVector(g)
            * self.population.fitnessVector(g) for g in genomes
            ]
        index_genome,rand = self.population.getWeightedRandom(
            rand,weights
            )
        index_other_genome,rand = self.population.getWeightedRandom(
            rand,weights
            )

        parents = [ genomes[index_genome], genomes[index_other_genome] ]
        recom_index = int( rand * self.population.num_loci )
        new_genome = parents[0][0:recom_index] + parents[1][recom_index:]
        if new_genome not in self.pathogens:
            self.acquirePathogen(new_genome)
