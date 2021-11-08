
"""Contains class Vector."""

import numpy as np
import copy as cp

class Vector(object):
    """Class defines vector entities to be infected by pathogens in model.

    These can infect hosts, the main entities in the model.

    Methods:
    copyState -- returns a slimmed-down version of the current vector state
    acquirePathogen -- adds given genome to this vector's pathogens
    infectHost -- infects given host with a sample of this vector's pathogens
    infectVector -- infects given vector with a sample of this vector's pathogens
    recover -- removes all infections
    die -- kills this vector
    birth -- add a new vector to population based on this vector
    applyTreatment -- removes all infections with genotypes susceptible to given
        treatment
    mutate -- mutate a single, random locus in a random pathogen
    recombine -- recombine two random pathogen genomes at random locus
    getWeightedRandomGenome -- returns index of element chosen from weights and
        given random number
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
            self.coefficient_index = population.coefficients_vectors.shape[0]

            population.coefficients_vectors = np.vstack( (
                population.coefficients_vectors,
                population.healthyCoefficientRow()
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
        sum_fitness_denom = self.sum_fitness if self.sum_fitness > 0 else 1
        self.population.coefficients_vectors[self.coefficient_index,:] = (
            self.population.coefficients_vectors[
                self.coefficient_index,:
                ]
            * old_sum_fitness / sum_fitness_denom ) + ( np.array([
                    # positions dependent on class constants
                0,
                self.population.contactVector(genome),
                self.population.receiveContactVector(genome),
                self.population.lethalityVector(genome),
                self.population.natalityVector(genome),
                self.population.recoveryVector(genome),
                self.population.migrationVector(genome),
                self.population.populationContactVector(genome),
                self.population.receivePopulationContactVector(genome),
                self.population.mutationVector(genome),
                self.population.recombinationVector(genome)
            ]) * self.pathogens[genome] / sum_fitness_denom )

        self.population.coefficients_vectors[
            self.coefficient_index,self.population.RECOMBINATION
            ] = self.population.coefficients_vectors[
                self.coefficient_index,self.population.RECOMBINATION
                ] * ( len(self.pathogens) > 1 )

        self.population.coefficients_vectors[
            self.coefficient_index,self.population.INFECTED
            ] = 1

        if self not in self.population.infected_vectors:
            self.population.infected_vectors.append(self)

    def infectHost(self, host):
        """Infect given host with a sample of this vector's pathogens.

        Each pathogen in the infector is sampled as present or absent in the
        inoculum by drawing from a Poisson distribution with a mean equal to the
        mean inoculum size of the organism being infected weighted by each
        genome's fitness as a fraction of the total in the infector as the
        probability of each trial (minimum 1 pathogen transfered). Each pathogen
        present in the inoculum will be added to the infected organism, if it
        does not have protection from the pathogen's genome. Fitnesses are
        computed for the pathogens' genomes in the infected organism, and the
        organism is included in the poplation's infected list if appropriate.

        Arguments:
        vector -- the vector to be infected (Vector)

        Returns:
        whether or not the model has changed state (Boolean)
        """

        changed = False

        genomes = list( self.pathogens.keys() )
        fitness_weights = [
            self.pathogens[g] / self.sum_fitness for g in genomes
            ]

        genomes_inoculated = np.unique( np.random.choice(
            genomes, p=fitness_weights,
            size=max(
                np.random.poisson( self.population.mean_inoculum_host ), 1
                )
            ) )
        for genome in genomes_inoculated:
            if genome not in host.pathogens.keys() and not any(
                    [ p in genome for p in host.protection_sequences ]
                    ):
                host.acquirePathogen(genome)
                changed = True

        return changed

    def infectVector(self, vector):
        """Infect given host with a sample of this vector's pathogens.

        Each pathogen in the infector is sampled as present or absent in the
        inoculum by drawing from a Poisson distribution with a mean equal to the
        mean inoculum size of the organism being infected weighted by each
        genome's fitness as a fraction of the total in the infector as the
        probability of each trial (minimum 1 pathogen transfered). Each pathogen
        present in the inoculum will be added to the infected organism, if it
        does not have protection from the pathogen's genome. Fitnesses are
        computed for the pathogens' genomes in the infected organism, and the
        organism is included in the poplation's infected list if appropriate.

        Arguments:
        vector -- the vector to be infected (Vector)

        Returns:
        whether or not the model has changed state (Boolean)
        """

        changed = False

        genomes = list( self.pathogens.keys() )
        fitness_weights = [
            self.pathogens[g] / self.sum_fitness for g in genomes
            ]

        genomes_inoculated = np.unique( np.random.choice(
            genomes, p=fitness_weights,
            size=max(
                np.random.poisson( self.population.mean_inoculum_vector ), 1
                )
            ) )
        for genome in genomes_inoculated:
            if genome not in vector.pathogens.keys() and not any(
                    [ p in genome for p in vector.protection_sequences ]
                    ):
                vector.acquirePathogen(genome)
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
                    seq = genome[
                        self.population.protection_upon_recovery_vector[0]
                        :self.population.protection_upon_recovery_vector[1]
                        ]
                    if seq not in self.protection_sequences:
                        self.protection_sequences.append(seq)

            self.pathogens = {}
            self.sum_fitness = 0
            self.population.coefficients_vectors[
                self.coefficient_index,:
                ] = self.population.healthyCoefficientRow()

            self.population.infected_vectors.remove(self)

    def die(self):
        """Add vector to population's dead list, remove it from alive ones."""

        if self in self.population.infected_vectors:
            self.population.infected_vectors.remove(self)

        for v in self.population.vectors[self.coefficient_index:]:
            v.coefficient_index -= 1

        self.population.coefficients_vectors = np.delete(
            self.population.coefficients_vectors, self.coefficient_index, 0
            )
        self.population.vectors.remove(self)

    def birth(self, rand):
        """Add vector to population based on this vector."""

        vector_list = self.population.addVectors(1)
        vector = vector_list[0]

        if self.population.vertical_transmission_vector > rand:
            self.infectVector(vector)

        if self.population.inherit_protection_vector > np.random.random():
            vector.protection_sequences = self.protection_sequences.copy()

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
                ] = self.population.healthyCoefficientRow()
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
        index_genome,rand = self.getWeightedRandomGenome( rand,weights )

        old_genome = genomes[index_genome]
        mut_index = int( rand * self.population.num_loci )
        if old_genome[mut_index] != self.population.CHROMOSOME_SEPARATOR:
            new_genome = old_genome[0:mut_index] + np.random.choice(
                list(self.population.possible_alleles[mut_index])
                ) + old_genome[mut_index+1:]
            if new_genome not in self.pathogens:
                self.acquirePathogen(new_genome)

            if new_genome not in self.population.model.global_trackers['genomes_seen']:
                self.population.model.global_trackers['genomes_seen'].append(new_genome)

    def recombine(self,rand):
        """Recombine two random pathogen genomes at random locus.

        Creates a new genotype from two random possible pathogens.
        """

        genomes = list( self.pathogens.keys() )
        weights = [
            self.population.recombinationVector(g)
            * self.population.fitnessVector(g) for g in genomes
            ]
        index_genome,rand = self.getWeightedRandomGenome( rand,weights )
        index_other_genome,rand = self.getWeightedRandomGenome( rand,weights )

        if index_genome != index_other_genome:
            num_evts = np.random.poisson( self.population.num_crossover_vector )
            loci = np.random.randint( 0, self.population.num_loci, num_evts )

            children = [ genomes[index_genome], genomes[index_other_genome] ]

            for l in loci:
                children[0] = children[0][0:l] + children[1][l:]
                children[1] = children[1][0:l] + children[0][l:]

            children = [
                genome.split(self.population.CHROMOSOME_SEPARATOR)
                for genome in children
                ]
            parent = np.random.randint( 0, 2, len( children[0] ) )

            children = [
                self.population.CHROMOSOME_SEPARATOR.join([
                    children[ parent[i] ][i]
                    for i in range( len( children[0] ) )
                    ]),
                self.population.CHROMOSOME_SEPARATOR.join([
                    children[ not parent[i] ][i]
                    for i in range( len( children[1] ) )
                    ])
                ]

            for new_genome in children:
                if new_genome not in self.pathogens:
                    self.acquirePathogen(new_genome)

                if new_genome not in self.population.model.global_trackers['genomes_seen']:
                    self.population.model.global_trackers['genomes_seen'].append(new_genome)

    def getWeightedRandomGenome(self, rand, r):
        """Returns index of element chosen from weights and given random number.

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
                return i, ( ( u - r_cum + e ) / e )
