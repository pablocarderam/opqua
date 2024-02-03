
"""Contains class Landscape."""

import copy as cp
import textdistance as td
import numpy as np
import pandas as pd

class Landscape(object):
    """Class defines a new fitness landscape mapping genomes and fitness.

    Methods:
    fitness -- evaluates fitness of genome using stored fitness function
    reduceGenome -- returns reduced genome with only fitness-relevant alleles
    establishmentRatePerPathogen -- returns establishment rate (per individual)
    makeNeighbor -- records one genome as having a second as a neighbor
    mapNeighbors -- recursively maps mutations in immediate neighborhood
    map -- maps and evaluates fitness of all relevant mutations
    save -- saves mutation network and fitness values
    load -- load mutation network and fitness values
    """

    def __init__(self, id=None, setup=None,
            fitnessFunc=None,
            population_threshold=None, selection_threshold=None,
            max_depth=None, allele_groups=None):
        """Create a new Landscape.

        Keyword arguments:
        id -- this landscape's indentifier (String)
        setup -- setup with associated parameters (Setup object, default None)
        fitnessFunc -- fitness function used to evaluate genomes (function
            taking a genome for argument and returning a fitness value >0,
            default None)
        population_threshold -- pathogen threshold under which drift is assumed
            to dominate (number >1, default None)
        selection_threshold -- selection coefficient threshold under which
            drift is assumed to dominate; related to population_threshold
            (number >1, default None)
        max_depth -- max number of mutations considered when evaluating
            establishment rates (integer >0, default None)
        allele_groups -- relevant alleles affecting fitness, each element
            contains a list of strings, each string contains a group of alleles
            that all have equivalent fitness behavior (list of lists of Strings)
        """
        super(Landscape, self).__init__()
        self.id = id
        self.setup = setup
        if setup is not None:
            setup.landscapes[id] = self

        self.fitnessFunc = setup.fitnessHost \
            if fitnessFunc is None else fitnessFunc
        self.population_threshold = setup.population_threshold \
            if population_threshold is None else population_threshold
        self.max_depth = setup.max_depth if max_depth is None else max_depth
        self.allele_groups = setup.allele_groups \
            if allele_groups is None else allele_groups
            # list where every item is a list of strings, each string containing
            # a group of alleles that all have equivalent fitness behavior

        if selection_threshold is not None:
            self.selection_threshold = selection_threshold
        elif population_threshold is not None:
            self.selection_threshold = 1 / population_threshold
        else:
            self.selection_threshold = setup.selection_threshold
            # selection threshold is defined
            # in terms of selection coefficient, i.e. selective advantage of
            # mutant over original genome with fitness=1

        self.equivalent_alleles = []
        self.num_equivalent_alleles = []
        self.relevant_loci = []
        self.allele_groups_reduced = []
        for i,locus in enumerate(setup.allele_groups):
            if len(locus) > 1:
                self.relevant_loci.append(i)
                reduced_alleles_locus = ''
                self.equivalent_alleles.append({})
                self.num_equivalent_alleles.append({})
                for group in locus:
                    reduced_alleles_locus = reduced_alleles_locus + group[0]
                    for allele in group:
                        self.equivalent_alleles[i][allele] = group[0]
                        self.num_equivalent_alleles[i][allele] = len(group)

                self.allele_groups_reduced.append(reduced_alleles_locus)

        self.fitness_values_reduced = {}
            # contain fitness values of reduced genomes

        self.mutation_network = {}
            # Contains all connections and rates:
            # Keys: reduced genomes. Values: dictionaries with
            #   'neighbors':list of neighboring reduced genomes
            #   'rates': list of corresponding establishment rates for neighbors
            #   'sum_rates': number with sum of all rates in previous list
        self.mapped = []

    def fitness(self,genome,genome_reduced=False, background_genome=None):
        """Fitness function used to evaluate genomes.

        Arguments:
        genome -- genome sequence to be evaluated (String)

        Keyword arguments:
        genome_reduced -- whether the genome given is reduced (Boolean)
        background_genome -- genome on which mutations are introduced, if
            sequence being evaluated is a reduced genome (String)

        Returns:
        Number >0 representing the fitness of genome g
        """
        if genome_reduced:
            if genome in self.fitness_values_reduced.keys():
                fitness = self.fitness_values_reduced[genome]
            elif background_genome is not None:
                full_genome = background_genome
                for i,locus in enumerate(self.relevant_loci):
                    full_genome = full_genome[0:locus] + genome[i] \
                        + full_genome[locus+1:]

                fitness = self.fitnessFunc(full_genome)
                self.fitness_values_reduced[genome] = fitness
            else:
                fitness = -1 # Error
        else:
            fitness = self.fitnessFunc(genome)

        return fitness

    def reduceGenome(self,genome):
        """Returns reduced genome with only fitness-relevant alleles.

        Arguments:
        genome -- full genome sequence to be evaluated (String)

        Returns:
        Reduced genome with only fitness-relevant alleles (String)
        """
        return ''.join( [
            self.equivalent_alleles[l][ genome[l] ] for l in self.relevant_loci
            ] )

    def establishmentRatePerPathogen(
            self, ancestor_fitness, mutant_fitness, distance, synonym_probs):
        """Returns establishment rate (calculated per individual pathogen).

        Arguments:
        ancestor_fitness -- fitness of ancestor genome (number)
        mutant_fitness -- fitness of mutant genome, >ancestor_fitness (number)
        distance -- Levenshtein distance between both genomes (integer)
        synonym_probs -- combinatorial probability of obtaining the exact
            mutations needed, given that same number of mutations happening
            (number)

        Returns:
        Rate at which the mutant fixates over the ancestor background (number)
        """
        return (
            np.power(
                self.setup.mutate_in_host / self.setup.generation_time, distance
                )
            * np.prod( synonym_probs )
            * ( ( mutant_fitness / ancestor_fitness ) - 1 )
                # this last factor is the selection coefficint
            )

    def makeNeighbor(
            self, ancestor_genome, ancestor_fitness,
            mutant_genome, mutant_fitness, distance, synonym_probs):
        """Records one genome as having a second as a neighbor

        Arguments:
        ancestor_genome -- ancestor reduced genome (String)
        ancestor_fitness -- fitness of ancestor genome (number)
        mutant_genome -- mutant reduced genome (String)
        mutant_fitness -- fitness of mutant genome (number)
        distance -- Levenshtein distance between both genomes (integer)
        synonym_probs -- combinatorial probability of obtaining the exact
            mutations needed, given that same number of mutations happening
            (number)
        """
        print('      Anc: '+ancestor_genome+', Mut: '+mutant_genome+', ',ancestor_fitness, mutant_fitness, distance, synonym_probs)
        if ancestor_genome not in self.mutation_network.keys():
            self.mutation_network[ ancestor_genome ] = {
                'neighbors':[], 'rates':[], 'sum_rates':0
                }

        establishment_rate = self.establishmentRatePerPathogen(
            ancestor_fitness, mutant_fitness, distance, np.prod(synonym_probs)
            )

        if ( mutant_genome
                not in self.mutation_network[ ancestor_genome ]['neighbors']):
            self.mutation_network[ ancestor_genome ]['neighbors'].append(
                mutant_genome
                )
            self.mutation_network[ ancestor_genome ]['rates'].append(
                establishment_rate
                )
            self.mutation_network[ ancestor_genome ]['sum_rates'] \
                += self.mutation_network[ ancestor_genome ]['rates'][-1]
        else:
            i = self.mutation_network[ ancestor_genome ]['neighbors'].index(
                mutant_genome
                )
            self.mutation_network[ ancestor_genome ]['sum_rates'] = (
                self.mutation_network[ ancestor_genome ]['sum_rates']
                - self.mutation_network[ ancestor_genome ]['rates'][i]
                + establishment_rate
                )
            self.mutation_network[
                ancestor_genome
                ]['rates'][i] = establishment_rate

    # def mapNeighbors(
    #         self, reduced_genome, fitness, background_genome, depth, ancestors):
    #     """Recursively maps and evaluates fitness of mutations in neighborhood.
    #
    #     Saves result in self.mutation_network property.
    #
    #     Recursive function!
    #
    #     Arguments:
    #     reduced_genome -- reduced genome used as parent for mutations (String)
    #     fitness -- fitness of background genome (number >1)
    #     background_genome -- full genome used as background (String)
    #     depth -- Levenshtein distance since last fitness increase (integer)
    #     ancestors -- List containing ancestor genomes and their associated
    #         properties; each element is a dictionary with:
    #             'genome': reduced ancestor genome (String)
    #             'fitness': fitness of reduced ancestor genome (String)
    #             'synonym_probs_fwd': combinatorial probability of obtaining the
    #                 exact mutations needed to get from ancestor to mutant, given
    #                 that same number of mutations happening (number)
    #             'synonym_probs_rev': combinatorial probability of obtaining the
    #                 exact mutations needed to get from mutant to ancestor, given
    #                 that same number of mutations happening (number)
    #     """
    #     print('Mapping '+reduced_genome)
    #     self.mapped.append(reduced_genome)
    #     for ancestor_distance,ancestor in enumerate(ancestors):
    #         selection_coefficient = ( ( fitness / ancestor['fitness'] ) - 1 )
    #             # of mutant vs. parent
    #
    #         if selection_coefficient > self.selection_threshold:
    #             depth = -1
    #             self.makeNeighbor(
    #                 ancestor['genome'], ancestor['fitness'],
    #                 reduced_genome, fitness,
    #                 ancestor_distance, ancestor['synonym_probs_fwd']
    #                 )
    #         elif ( -1 * selection_coefficient > self.selection_threshold ):
    #             self.makeNeighbor(
    #                 reduced_genome, fitness,
    #                 ancestor['genome'], ancestor['fitness'],
    #                 ancestor_distance, ancestor['synonym_probs_rev']
    #                 )
    #
    #     if depth < self.max_depth:
    #         depth += 1
    #         new_ancestors = [{
    #                 'genome':reduced_genome,'fitness':fitness,
    #                 'synonym_probs_fwd':1,'synonym_probs_rev':1
    #             }] + cp.deepcopy(ancestors[0:-2])
    #         for locus,locus_alleles in enumerate( self.allele_groups_reduced ):
    #             for allele in locus_alleles:
    #                 if allele != reduced_genome[locus]:
    #                     mutant = reduced_genome[0:locus] \
    #                         + allele + reduced_genome[locus+1:]
    #                     print('Scanning mutant: '+mutant)
    #                     if mutant not in self.mapped:
    #                         mutant_fitness = self.fitness(
    #                             mutant, genome_reduced=True,
    #                             background_genome=background_genome
    #                             )
    #                         for ancestor in new_ancestors:
    #                             ancestor['synonym_probs_fwd'] = (
    #                                 ancestor['synonym_probs_fwd']
    #                                 * self.num_equivalent_alleles[locus][
    #                                     allele
    #                                     ] / len( self.allele_groups[locus] )
    #                                 )
    #                             ancestor['synonym_probs_rev'] = (
    #                                 ancestor['synonym_probs_rev']
    #                                 * self.num_equivalent_alleles[locus][
    #                                     reduced_genome[locus]
    #                                     ] / len( self.allele_groups[locus] )
    #                                 )
    #                         self.mapNeighbors(
    #                             mutant, mutant_fitness, background_genome,
    #                             depth, new_ancestors
    #                             )
    def evaluateNeighbors(
            self, reduced_genome, fitness, background_genome, depth,
            synonym_probs_fwd, synonym_probs_rev, ancestor, ancestor_fitness,
            loci_mutated):
        """Recursively evaluates fitness of mutations in neighborhood.

        Saves result in self.mutation_network property.

        Recursive function!

        Arguments:
        reduced_genome -- reduced genome used as parent for mutations (String)
        fitness -- fitness of background genome (number >1)
        background_genome -- full genome used as background (String)
        depth -- Levenshtein distance since last fitness increase (integer)
        ancestors -- List containing ancestor genomes and their associated
            properties; each element is a dictionary with:
                'genome': reduced ancestor genome (String)
                'fitness': fitness of reduced ancestor genome (String)
                'synonym_probs_fwd': combinatorial probability of obtaining the
                    exact mutations needed to get from ancestor to mutant, given
                    that same number of mutations happening (number)
                'synonym_probs_rev': combinatorial probability of obtaining the
                    exact mutations needed to get from mutant to ancestor, given
                    that same number of mutations happening (number)
        """
        print('   Scanning from: '+reduced_genome+', depth: '+str(depth),', ancestor: '+ancestor+', loci mutated: '+str(loci_mutated))
        first_neighbors = {}
        for locus,locus_alleles in enumerate( self.allele_groups_reduced ):
            if locus not in loci_mutated:
                for allele in locus_alleles:
                    if allele != reduced_genome[locus]:
                        mutant = reduced_genome[0:locus] \
                            + allele + reduced_genome[locus+1:]
                        print('      Scanning mutant: '+mutant)

                        if mutant != ancestor:
                            mutant_fitness = self.fitness(
                                mutant, genome_reduced=True,
                                background_genome=background_genome
                                )
                            selection_coefficient = (
                                ( mutant_fitness / ancestor_fitness ) - 1
                                )
                                # of mutant vs. ancestor

                            if depth == 1:
                                first_neighbors[mutant] = {
                                    'fitness':mutant_fitness,
                                    'more_fit':False
                                    }

                            synonym_probs_fwd_new = synonym_probs_fwd + [
                                self.num_equivalent_alleles[locus][
                                    allele
                                    ] / len( self.allele_groups[locus] )
                                ]
                            synonym_probs_rev_new = synonym_probs_rev + [
                                self.num_equivalent_alleles[locus][
                                    reduced_genome[locus]
                                    ] / len( self.allele_groups[locus] )
                                ]

                            if selection_coefficient > self.selection_threshold:
                                self.makeNeighbor(
                                    ancestor, ancestor_fitness,
                                    mutant, mutant_fitness,
                                    depth,synonym_probs_fwd_new
                                    )
                                if depth == 1:
                                    first_neighbors[mutant]['more_fit'] = True

                            elif ( -1 * selection_coefficient
                                    > self.selection_threshold ):
                                self.makeNeighbor(
                                    mutant, mutant_fitness,
                                    ancestor, ancestor_fitness,
                                    depth,synonym_probs_rev_new
                                    )

                            if depth < self.max_depth:
                                self.evaluateNeighbors(
                                    mutant, mutant_fitness, background_genome,
                                    depth+1,
                                    synonym_probs_fwd_new, synonym_probs_rev_new,
                                    ancestor, ancestor_fitness,
                                    loci_mutated + [locus]
                                    )

        return first_neighbors

    def mapNeighbors(
            self, reduced_genome, fitness, background_genome, depth):
        """Recursively maps and evaluates fitness of mutations in neighborhood.

        Saves result in self.mutation_network property.

        Recursive function!

        Arguments:
        reduced_genome -- reduced genome used as parent for mutations (String)
        fitness -- fitness of background genome (number >1)
        background_genome -- full genome used as background (String)
        depth -- Levenshtein distance since last fitness increase (integer)
        ancestors -- List containing ancestor genomes and their associated
            properties; each element is a dictionary with:
                'genome': reduced ancestor genome (String)
                'fitness': fitness of reduced ancestor genome (String)
                'synonym_probs_fwd': combinatorial probability of obtaining the
                    exact mutations needed to get from ancestor to mutant, given
                    that same number of mutations happening (number)
                'synonym_probs_rev': combinatorial probability of obtaining the
                    exact mutations needed to get from mutant to ancestor, given
                    that same number of mutations happening (number)
        """
        print('Mapping '+reduced_genome+', depth: '+str(depth))
        self.mapped.append(reduced_genome)
        depth += 1
        first_neighbors = self.evaluateNeighbors(
            reduced_genome, fitness, background_genome, depth, [1], [1],
            reduced_genome, fitness, []
            )
        print('First neighbors of '+reduced_genome+': '+','.join(first_neighbors.keys()))
        for mutant in first_neighbors.keys():
            if mutant not in self.mapped:
                if first_neighbors[mutant]['more_fit']:
                    depth = 0

                if depth < self.max_depth:
                    self.mapNeighbors(
                        mutant, first_neighbors[ mutant ]['fitness'],
                        background_genome, depth
                        )

    def map(self,seed_genomes):
        """Maps and evaluates relevant mutations given object parameters

        Saves result in self.mutation_network property.

        Arguments:
        seed_genomes -- genome or list of genomes used as background for
            mutations (String or list of Strings)
        """
        print('Landscape mapping starting.')

        if not isinstance(seed_genomes,list):
            seed_genomes = [seed_genomes]

        for seed_genome in seed_genomes:
            fitness = self.fitness(seed_genome)
            reduced_genome = self.reduceGenome(seed_genome)
            self.mapped = []

            self.mapNeighbors( reduced_genome, fitness, seed_genome, 0 )

        print('Landscape mapping complete.')

    def save(save_to_file):
        """Saves mutation network and fitness values stored in mutation_network

        CSV format has the following columns:
        Genome: reduced genome
        Neighbors: list of neighboring reduced genomes, separated by semicolons
        Rates: list of corresponding establishment rates for neighbors,
            separated by semicolons
        Sum_rates: number with sum of all rates in previous list

        Arguments:
        save_to_file -- file path and name to save model data under (String)
        """
        print('Saving landscape to file...')

        out = 'Genome,Neighbors,Rates,Sum_rates\n'
        for genome in self.mutation_network:
            out = out + genome + ',' + ';'.join(
                self.mutation_network[genome]['neighbors']
                ) + ',' + ';'.join(
                self.mutation_network[genome]['rates']
                ) + ',' + str( self.mutation_network[genome]['sum_rates'] )+'\n'

        file = open(save_to_file,'w')
        file.write(out)
        file.close()

        print('...file saved.')

    def load(file):
        """Loads mutation network and fitness from file path

        CSV format has the following columns:
        Genome: reduced genome
        Neighbors: list of neighboring reduced genomes, separated by semicolons
        Rates: list of corresponding establishment rates for neighbors,
            separated by semicolons
        Sum_rates: number with sum of all rates in previous list

        Arguments:
        file -- file path and name to save model data under (String)
        """
        print('Loading landscape from file...')

        df = pd.load_csv(file)

        self.mutation_network = {}
        for i,row in df.iterrows():
            self.mutation_network[ row['Genome'] ] : {
                'neighbors' : row['Neighbors'].split(';'),
                'rates' : [ float(r) for r in row['Rates'].split(';') ],
                'sum_rates' : float( row['Sum_rates'] )
                }

        print('...landscape loaded.')
