
"""Contains class Landscape."""

import math as ma
import copy as cp
import textdistance as td
import numpy as np
import pandas as pd
import scipy.special as sp_spe
import scipy.stats as sp_sta
import scipy.integrate as sp_int
import itertools as itt

class Landscape(object):
    """Class defines a new fitness landscape mapping genomes and fitness.

    Methods:
    fitness -- evaluates fitness of genome using stored fitness function
    reduceGenome -- returns reduced genome with only fitness-relevant alleles
    survivalProbabilities -- returns array of survival probabilities per
        generation for mutant
    establishmentRates -- returns establishment rate (per individual)
    makeNeighbor -- records one genome as having a second as a neighbor
    evaluateNeighbors -- recursively evaluates fitness of neighboring mutations
    mapNeighbors -- recursively maps mutations in immediate neighborhood
    map -- maps and evaluates fitness of all relevant mutations
    save -- saves mutation network and fitness values
    load -- load mutation network and fitness values
    """

    def __init__(self, id=None, setup=None,
            fitnessFunc=None, mutationFunc=None, mutate=None,
            generation_time=None, population_threshold=None,
            selection_threshold=None, max_depth=None, allele_groups=None,
            population_size=None, max_generations_survival=None):
        """Create a new Landscape.

        Keyword arguments:
        id -- this landscape's indentifier (String)
        setup -- setup with associated parameters (Setup object, default None)
        fitnessFunc -- fitness function used to evaluate genomes (function
            taking a genome for argument and returning a fitness value >0,
            default None)
        mutationFunc -- mutation function used to evaluate genomes (function
            taking a genome for argument and returning a mutation coefficient>0,
            default None)
        mutate -- mutation rate per generation (number>0, default None)
        generation_time -- time between pathogen generations (number>0, default
            None)
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
        population_size -- intrahost population (integer)
        max_generations_survival -- number of generations used to compute rates
            and probabilities of genotype emergence (integer)
        """
        super(Landscape, self).__init__()
        self.id = id
        self.setup = setup
        if setup is not None:
            setup.landscapes[id] = self

        self.fitnessFunc = setup.fitnessHost \
            if fitnessFunc is None else fitnessFunc
        self.mutationFunc = setup.mutationHost \
            if mutationFunc is None else mutationFunc
        self.mutate = setup.mutate_in_host if mutate is None else mutate
        self.generation_time = setup.generation_time_host \
            if generation_time is None else generation_time
        self.population_threshold = setup.population_threshold_host \
            if population_threshold is None else population_threshold
        self.population_size = setup.steady_pathogen_population_host \
            if population_size is None else population_size
        self.max_generations_survival = int(setup.max_generations_survival_host) \
            if max_generations_survival is None else int(max_generations_survival)
        self.max_depth = setup.max_depth_host if max_depth is None else max_depth
        self.allele_groups = setup.allele_groups_host \
            if allele_groups is None else allele_groups
            # list where every item is a list of strings, each string containing
            # a group of alleles that all have equivalent fitness behavior

        if selection_threshold is not None:
            self.selection_threshold_host = selection_threshold
        elif population_threshold is not None:
            self.selection_threshold = 1 / population_threshold
        else:
            self.selection_threshold = setup.selection_threshold_host
            # selection threshold is defined
            # in terms of selection coefficient, i.e. selective advantage of
            # mutant over original genome with fitness=1

        self.equivalent_alleles = []
        self.num_equivalent_alleles = []
        self.relevant_loci = []
        self.allele_groups_reduced = []
        self.total_alleles = []
        for i,locus in enumerate(self.allele_groups):
            if len(locus) > 1:
                self.relevant_loci.append(i)
                reduced_alleles_locus = ''
                self.equivalent_alleles.append({})
                self.num_equivalent_alleles.append({})
                self.total_alleles.append(0)
                for group in locus:
                    reduced_alleles_locus = reduced_alleles_locus + group[0]
                    for allele in group:
                        self.equivalent_alleles[i][allele] = group[0]
                        self.num_equivalent_alleles[i][allele] = len(group)
                        self.total_alleles[i] += 1

                self.allele_groups_reduced.append(reduced_alleles_locus)

        self.fitness_values_reduced = {}
            # contain fitness values of reduced genomes
        self.mutation_rates_reduced = {}
            # contain mutation rates of reduced genomes

        self.mutation_network = {}
            # Contains all connections and rates:
            # Keys: reduced genomes. Values: dictionaries with
            #   'neighbors':list of neighboring reduced genomes
            #   'rates': list of lists, each inner list has the list of
            #       establishment rates for the corresponding neighbor in
            #       'neighbors', each rate in inner list corresponds to rate at
            #       the mutant cloud size corresponding to the inner index
            #   'sum_rates': number with sum of all rates in previous list
        self.mapped = []

        # Fixed evolutionary parameters
        self.max_cloud_size_factor = 2
            # how many times bigger than max depth we will consider
        self.alpha = 4 # for Gamma distribution of fixation waiting times
        self.scaling_coefficient = 0.67
            # for Gamma distribution of fixation waiting times
        self.m_param = -1.12 # for correction of mean fixation waiting time
        self.k_param = 1.27 # for correction of mean fixation waiting time

        # Precalculate PMF values for use below
        self.poisson_pmf_vals = np.zeros(
            ( self.max_depth, self.max_depth * self.max_cloud_size_factor )
            )
        for distance in range( self.max_depth ):
            for d_cloud_em in range(self.max_depth *self.max_cloud_size_factor):
                # poisson_pmf_vals[d_cloud_em] = sp_sta.poisson.pmf(
                #     distance, d_cloud_em
                #     )
                self.poisson_pmf_vals[ distance, d_cloud_em ] = (
                    np.power( d_cloud_em, distance ) * np.exp( -d_cloud_em )
                    / ma.factorial( distance )
                    )

        # Will store precalculated PMF values for use below
        self.seed_mutation_rate = 0
        self.geom_pmf_vals_seed = []

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

    def mutationRate(self,genome,genome_reduced=False, background_genome=None):
        """Mutation function used to evaluate genomes.

        Arguments:
        genome -- genome sequence to be evaluated (String)

        Keyword arguments:
        genome_reduced -- whether the genome given is reduced (Boolean)
        background_genome -- genome on which mutations are introduced, if
            sequence being evaluated is a reduced genome (String)

        Returns:
        Number >0 representing the mutation rate (per genome) of genome g
        """
        if genome_reduced:
            if genome in self.mutation_rates_reduced.keys():
                mutation_rate = self.mutation_rates_reduced[genome]
            elif background_genome is not None:
                full_genome = background_genome
                for i,locus in enumerate(self.relevant_loci):
                    full_genome = full_genome[0:locus] + genome[i] \
                        + full_genome[locus+1:]

                mutation_rate = (
                    self.mutationFunc(full_genome)
                    * self.mutate * len( self.relevant_loci )
                    )
                self.mutation_rates_reduced[genome] = mutation_rate
            else:
                mutation_rate = -1 # Error
        else:
            mutation_rate = (
                self.mutationFunc(genome)
                * self.mutate * len( self.relevant_loci )
                )

        return mutation_rate

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

    def survivalProbabilities(
            self, selection_coefficient, max_generations_survival,
            population_size, population_size_threshold=1e4):
        """Returns array of survival probabilities per generation for mutant.

        Arguments:
        selection_coefficient -- between mutant and background (number)
        max_generations_survival -- max generations computed (integer)
        population_size -- population size considered (integer)
        population_size_threshold -- population over which approximation is used
            (integer)

        Returns:
        array of survival probabilities per generation for mutant (array)
        """
        probs = np.zeros( max_generations_survival + 1 )
        if population_size > population_size_threshold:
            for i in range( 1, max_generations_survival + 1 ):
                probs[i] = np.exp(
                    ( 1 + selection_coefficient ) * ( probs[i-1] - 1 )
                    ) # Fisher (unfortunately) 1922

        else:
            for i in range( 1, max_generations_survival + 1 ):
                probs[i] = np.power(
                    1 - (
                        ( 1 + selection_coefficient ) * ( 1 - probs[i-1] )
                        / population_size
                        ),
                    population_size
                    )

        return 1-probs[1:]

    def establishmentRates(
            self, ancestor_fitness, mutant_fitness, distance,
            mutations, ancestor_genome, background_genome):
        """Returns establishment rates (calculated per individual pathogen).

        One establishment rate per mutant cloud distance is calculated.

        Arguments:
        ancestor_fitness -- fitness of ancestor genome (number)
        mutant_fitness -- fitness of mutant genome, >ancestor_fitness (number)
        distance -- Hamming distance between both genomes (integer)
        mutations -- list of mutations  between the ancestor genome and mutant,
            in format "WT allele"+"position"+"mutant allele" (e.g. A334K) (list)
        ancestor_genome -- reduced ancestor genome sequence (String)
        background_genome -- full background genome that reduced ancestor was
            extracted from (String)

        Returns:
        Array with rates at which the mutant fixates over the ancestor
        background for each mutant cloud radius (corresponding to element index)
        (array of numbers)
        """
        # print('            Establishment rates: '+str([ancestor_genome,mutations]))

        ancestor_mutation_rate = self.mutationRate( ancestor_genome )
        if ancestor_mutation_rate == self.seed_mutation_rate:
            geom_pmf_vals_anc = self.geom_pmf_vals_seed
        else:
            max_generations_until_emergence = int( np.ceil(
                self.max_depth * self.max_cloud_size_factor
                / ancestor_mutation_rate
                ) )
            geom_pmf_vals_anc = np.zeros(max_generations_until_emergence)
            for x in range(max_generations_until_emergence):
                # geom_pmf_vals_anc[x] = sp_sta.geom.pmf( x+1,ancestor_mutation_rate )
                geom_pmf_vals_anc[x] = np.power(
                    1-ancestor_mutation_rate, x
                    ) * x

        rates = [ np.zeros(1) ] * ( self.max_depth * self.max_cloud_size_factor )

        q_em = 0 # population fraction with this mutant genotype (emergence)
        p_fix = 0 # probability of this mutant fixating (establishing at 1:1)
            # once it has emerged
        q_fix = 0 # fraction of fixation events occurring during this current
            # time frame (before mutant cloud distance is advanced)
        t_len = 0 # length of current time frame (before mutant cloud distance
            # is advanced)
        selection_coefficient = (
            ( mutant_fitness / ancestor_fitness ) - 1
            )

        t_len = self.generation_time / ancestor_mutation_rate

        p_fix = (
            ( - np.exp( -2 * selection_coefficient ) + 1 )
            / ( - np.exp(
                -2 * self.population_size * selection_coefficient
                ) + 1 )
            ) # Kimura diffusion approximation (1957, 1962)

        p_cm = (
            np.power( 1/len( self.relevant_loci ), distance )
            * np.prod([
                self.num_equivalent_alleles[ int(mut[1:-1]) ][ mut[-1] ]
                / self.total_alleles[ int(mut[1:-1]) ]
                for mut in mutations
                ])
            )
            # probability of getting correct combination of mutations
            # given correct number of mutations

        mut_permutations = list( itt.permutations( mutations ) )

        mean_fix_t = (
                # mean fixation (establishment to 50%) time of mutant
            np.log( self.population_size )
            / np.log(selection_coefficient+1)
                # this is from Tami Liebermann's class, I think derives
                # from Kimura 1957, 1962 again
            + np.exp(
                self.m_param * np.log( selection_coefficient ) + self.k_param
                )
                # this is a correction factor that matches up, derived
                # from simulation, parameter weights fitted
            )
        beta = selection_coefficient / self.scaling_coefficient
            # beta to be used in Gamma/Erland PDF; these values fitted
            # from simulation

        for d_cloud_current in range( self.max_depth * self.max_cloud_size_factor ):
                # d_cloud_current describes mutant cloud size in units of
                # mutations
            rates[d_cloud_current] = np.zeros(d_cloud_current)
            # print('              Cloud current: '+str(d_cloud_current))
            for d_cloud_em in range( d_cloud_current ):
                # print('                Cloud emergence: '+str(d_cloud_em))
                    # d_cloud_em describes mutant cloud size in units of
                    # mutations at moment in which mutation arose
                q_nm = 0 # population fraction with correct number of mutations
                if d_cloud_em == 0:
                        # if mutant cloud just started, assume variation from 1
                        # replication cycle
                    q_nm = np.power( ancestor_mutation_rate, distance )
                else:
                        # otherwise, variation is set by Poisson PMF mutant
                        # cloud with mean=cloud size and k=distance of mutant
                        # of interest
                    # print(d_cloud_em)
                    q_nm = self.poisson_pmf_vals[ distance, d_cloud_em ]

                sum_w = 0
                    # sum of likelihood weights for every available mutation
                    # path leading to mutant

                weight_emergence_permutations =np.zeros( len(mut_permutations) )
                generations_until_emergence = int( np.ceil(
                    ( d_cloud_current - d_cloud_em )
                    / ancestor_mutation_rate
                    ) )
                mutation_probabilities_anc = np.array([
                    geom_pmf_vals_anc[x]
                    for x in range(1, generations_until_emergence+1)
                    ])
                survival_probabilities_anc \
                        = self.survivalProbabilities(
                    1, generations_until_emergence, self.population_size
                    )

                for i,path in enumerate(mut_permutations):
                    # print('Rates for path '+str(path))
                    prev_intermediate = ancestor_genome
                    weight_intermediates = np.zeros( len(path) )
                    for j,mutation in enumerate(path):
                        intermediate = prev_intermediate[
                            0:int( mutation[1:-1] )
                            ] + mutation[-1] \
                            + prev_intermediate[ (int( mutation[1:-1] ) + 1): ]
                        # print('                  Intermediate: '+str(intermediate))
                        intermediate_fitness = self.fitness(
                            intermediate, genome_reduced=True,
                            background_genome=background_genome
                            )
                        intermediate_selection_coefficient = (
                            ( intermediate_fitness / ancestor_fitness ) - 1
                            ) # of intermediate vs. ancestor
                        if intermediate_selection_coefficient <= 0:
                                # if not detrimental, then we assume mutant is
                                # governed by drift (because this is prior to
                                # establishment)
                            intermediate_mutation_rate = self.mutationRate(
                                intermediate, genome_reduced=True,
                                background_genome=background_genome
                                )
                            if ( intermediate_mutation_rate
                                    == ancestor_mutation_rate ):
                                mutation_probabilities_mut \
                                    = mutation_probabilities_anc
                            else:
                                mutation_probabilities_mut = np.array([
                                    # sp_sta.geom.pmf(
                                    #     x+1,intermediate_mutation_rate
                                    #     )
                                    np.power( 1-intermediate_mutation_rate,x )*x
                                        # change this to save mut probabilities
                                        # in class pmf data structure and
                                        # reaccess here if possible? perhaps not
                                    for x in range(
                                        1, generations_until_emergence+1
                                        )
                                    ])
                            survival_probabilities_mut \
                                    = self.survivalProbabilities(
                                intermediate_selection_coefficient,
                                generations_until_emergence,
                                self.population_size
                                )

                            weight_intermediates[j] = (
                                np.sum( np.multiply(
                                    mutation_probabilities_mut,
                                    survival_probabilities_mut
                                    ) )
                                / np.sum( np.multiply(
                                    mutation_probabilities_anc,
                                    survival_probabilities_anc
                                    ) )
                                )
                        else:
                            weight_intermediates[j] = 1

                        prev_intermediate = intermediate

                    weight_emergence_permutations[i] \
                        = np.prod(weight_intermediates)

                    # print('Weight for path: '+str(weight_emergence_permutations[i]))
                sum_w = weight_emergence_permutations.sum()

                q_em = q_nm * p_cm * sum_w
                    # fraction that have emerged depend on fraction with correct
                    # number of mutations, probability of having correct
                    # combination of mutations, and weighted likelihood of
                    # mutant lineages surviving and acquiring those mutations
                    # (weighted relative to neutral drift)

                q_fix = sp_int.quad( # numerically computes a definite integral
                    lambda x: max( 0,
                        np.power( x-mean_fix_t, self.alpha-1 )
                        * np.exp( -beta * ( x-mean_fix_t ) )
                        * np.power( beta, self.alpha )
                        / ma.gamma( self.alpha )
                        ),
                        # integrating Gamma PDF describing distribution of
                        # mutant fixation wait times
                    ( d_cloud_current-d_cloud_em-1 ) / (ancestor_mutation_rate),
                    ( d_cloud_current-d_cloud_em ) / ( ancestor_mutation_rate )
                    # time interval boundaries of integration given by times of
                    # current mutant cloud size being considered
                    )[0] # first argument is integration, second is error

                rates[d_cloud_current][d_cloud_em] = (
                    q_em * p_fix * q_fix / t_len
                    )

            rates[d_cloud_current] = rates[d_cloud_current].sum()

        return np.array(rates)

    def makeNeighbor(
            self, ancestor_genome, ancestor_fitness,
            mutant_genome, mutant_fitness, distance, mutations,
            background_genome):
        """Records one genome as having a second as a neighbor

        Arguments:
        ancestor_genome -- ancestor reduced genome (String)
        ancestor_fitness -- fitness of ancestor genome (number)
        mutant_genome -- mutant reduced genome (String)
        mutant_fitness -- fitness of mutant genome (number)
        distance -- Hamming distance between both genomes (integer)
        mutations -- list of mutations  between the ancestor genome and mutant,
            in format "WT allele"+"position"+"mutant allele" (e.g. A334K) (list)
        background_genome -- full background genome that reduced ancestor was
            extracted from (String)
        """
        # print('          Making neighbor: '+str([ancestor_genome,mutant_genome]))
        # print('      Anc: '+ancestor_genome+', Mut: '+mutant_genome+', ',ancestor_fitness, mutant_fitness, distance, synonym_probs)
        if ancestor_genome not in self.mutation_network.keys():
            self.mutation_network[ ancestor_genome ] = {
                'neighbors':[], 'rates':[], 'sum_rates':0,
                'fitness':ancestor_fitness
                }
        if mutant_genome not in self.mutation_network.keys():
            self.mutation_network[ mutant_genome ] = {
                'neighbors':[], 'rates':[], 'sum_rates':0,
                'fitness':mutant_fitness
                }

        # print('Establishment rates from '+ancestor_genome+' to '+mutant_genome+' (d='+str(distance)+', s='+str((mutant_fitness/ancestor_fitness)-1)+')')
        establishment_rates = self.establishmentRates(
            ancestor_fitness, mutant_fitness, distance, mutations,
            ancestor_genome, background_genome
            )

        if ( mutant_genome
                not in self.mutation_network[ ancestor_genome ]['neighbors']):
            self.mutation_network[ ancestor_genome ]['neighbors'].append(
                mutant_genome
                )
            self.mutation_network[ ancestor_genome ]['rates'].append(
                establishment_rates
                )
            self.mutation_network[ ancestor_genome ]['sum_rates'] \
                += self.mutation_network[ ancestor_genome ]['rates'][-1].sum()
        else:
            i = self.mutation_network[ ancestor_genome ]['neighbors'].index(
                mutant_genome
                )
            self.mutation_network[ ancestor_genome ]['sum_rates'] = (
                self.mutation_network[ ancestor_genome ]['sum_rates']
                - self.mutation_network[ ancestor_genome ]['rates'][i].sum()
                + establishment_rates.sum()
                )
            self.mutation_network[
                ancestor_genome
                ]['rates'][i] = establishment_rates

    def evaluateNeighbors(
            self, reduced_genome, fitness, background_genome, depth,
            mutations, ancestor, ancestor_fitness, loci_mutated):
        """Recursively evaluates fitness of mutations in neighborhood.

        Saves result in self.mutation_network property.

        Recursive function!

        Arguments:
        reduced_genome -- reduced genome used as parent for mutations (String)
        fitness -- fitness of background genome (number >1)
        background_genome -- full genome used as background (String)
        depth -- Hamming distance since last fitness increase (integer)
        mutations -- list of mutations  between the ancestor genome and mutant,
            in format "WT allele"+"position"+"mutant allele" (e.g. A334K) (list)
        ancestor -- ancestor genome from which neighbors are evaluated (String)
        ancestor_fitness -- fitness of ancestor (number)
        loci_mutated -- list of indexes of positions already mutated in this
            evaluation run (List of integers)

        Returns:
        dictionary with keys=all first neighbor genomes, values=dictionary with
        'fitness':fitness of neighbor, 'more_fit':whether or not this mutant is
        more fit than the ancestor being evaluated (above the selection
        threshold).
        """
        # print('   Scanning from: '+reduced_genome+', depth: '+str(depth),', ancestor: '+ancestor+', loci mutated: '+str(loci_mutated))
        depth += 1
        first_neighbors = {}
        # print('    Evaluating: '+reduced_genome)
        for locus,locus_alleles in enumerate( self.allele_groups_reduced ):
            # print('      Locus: '+str([locus,locus_alleles]))
            if locus not in loci_mutated:
                for allele in locus_alleles:
                    # print('        Allele: '+str(allele))
                    if allele != reduced_genome[locus]:
                        mutant = reduced_genome[0:locus] \
                            + allele + reduced_genome[locus+1:]
                        # print('      Scanning mutant: '+mutant)

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

                            if selection_coefficient > self.selection_threshold:
                                self.makeNeighbor(
                                    ancestor, ancestor_fitness,
                                    mutant, mutant_fitness,
                                    depth, mutations+[
                                        reduced_genome[locus] + str(locus)
                                        + allele
                                        ],
                                    background_genome
                                    )
                                if depth == 1:
                                    first_neighbors[mutant]['more_fit'] = True

                            elif ( -1 * selection_coefficient
                                    > self.selection_threshold ):
                                flipped_mutations = [
                                    mutation[-1] + mutation[1:-1] + mutation[0]
                                    for mutation in mutations
                                    ]
                                self.makeNeighbor(
                                    mutant, mutant_fitness,
                                    ancestor, ancestor_fitness,
                                    depth, flipped_mutations+[
                                        allele + str(locus)
                                        + reduced_genome[locus]
                                        ],
                                    background_genome
                                    )

                            if depth < self.max_depth:
                                self.evaluateNeighbors(
                                    mutant, mutant_fitness, background_genome,
                                    depth,
                                    mutations+[
                                        reduced_genome[locus] + str(locus)
                                        + allele
                                        ],
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
        depth -- Hamming distance since last fitness increase (integer)
        """
        # print('Mapping '+reduced_genome+', depth: '+str(depth))
        self.mapped.append(reduced_genome)
        depth += 1
        print('  Evaluating: '+reduced_genome)
        first_neighbors = self.evaluateNeighbors(
            reduced_genome, fitness, background_genome, 0, [],
            reduced_genome, fitness, []
            )
        # print('First neighbors of '+reduced_genome+': '+','.join(first_neighbors.keys()))
        for mutant in first_neighbors.keys():
            if mutant not in self.mapped:
                if first_neighbors[mutant]['more_fit']:
                    depth = 0

                if depth < self.max_depth:
                    print('  Mapping: '+mutant)
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

            self.seed_mutation_rate = self.mutationRate( seed_genome )
            original_max_generations_until_emergence = int( np.ceil(
                self.max_depth * self.max_cloud_size_factor
                / self.seed_mutation_rate
                ) )
            self.geom_pmf_vals_seed = np.zeros(
                original_max_generations_until_emergence
                )
            for x in range(original_max_generations_until_emergence):
                self.geom_pmf_vals_seed[x] = np.power(
                    1-self.seed_mutation_rate, x
                    ) * x

            print('Reduced seed: '+seed_genome)
            self.mapNeighbors( reduced_genome, fitness, seed_genome, 0 )

        print('Landscape mapping complete.')

    def save(self,save_to_file):
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

        out = 'Genome,Neighbors,Rates,Sum_rates,Fitness\n'
        for genome in self.mutation_network:
            out = out + genome + ',' + ';'.join(
                self.mutation_network[genome]['neighbors']
                ) + ',' + ';'.join(
                [
                    '_'.join( [ str(r) for r in rates ] )
                    for rates in self.mutation_network[genome]['rates']
                    ]
                ) + ',' + str( self.mutation_network[genome]['sum_rates'] ) \
                + ',' + str( self.mutation_network[genome]['fitness'] )+'\n'

        file = open(save_to_file,'w')
        file.write(out)
        file.close()

        print('...file saved.')

    def load(self,file):
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

        df = pd.read_csv(file).fillna('')

        self.mutation_network = {}
        for i,row in df.iterrows():
            self.mutation_network[ row['Genome'] ] = {
                'neighbors' : row['Neighbors'].split(';'),
                'rates' : [
                    [ r for r in rates.split('_') ]
                    for rates in row['Rates'].split(';')
                    ],
                'sum_rates' : float( row['Sum_rates'] ),
                'fitness' : float( row['Fitness'] )
                }
            if len( self.mutation_network[ row['Genome'] ]['neighbors'][0] ) >0:
                self.mutation_network[ row['Genome'] ]['rates'] \
                     = np.asarray(
                     self.mutation_network[ row['Genome'] ]['rates']
                     ).astype(float)
            else:
                self.mutation_network[ row['Genome'] ]['neighbors'] = []
                self.mutation_network[ row['Genome'] ]['rates'] = []

        print('...landscape loaded.')
