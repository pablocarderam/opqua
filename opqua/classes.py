
import numpy as np
import copy as cp

class Host(object):
    """docstring for Host."""

    def __init__(self, population, id):
        super(Host, self).__init__()
        self.id = id
        self.pathogens = {}
        self.population = population
        self.sum_fitness = 0
        self.protection_sequences = []


    def infectVector(self, vector):
        """ Infects given vector """

        changed = False
        for genome,fitness in self.pathogens.items():
            if genome not in vector.pathogens.keys() and genome not in vector.protection_sequences and np.random.binomial( self.population.inoculum_host, self.population.inoculation_rate_host * fitness / self.sum_fitness, 1 ) > 0:
                vector.pathogens[genome] = vector.population.fitnessVector(genome)
                vector.sum_fitness += vector.pathogens[genome]
                vector.population.num_infected_vectors += 1
                changed = True

            rand = np.random.random()
            if self.population.mutate_in_vector > rand and changed:
                self.population.mutate(vector)
                rand = rand / self.population.mutate_in_vector
            else:
                rand = ( rand - self.population.mutate_in_vector ) / ( 1 - self.population.mutate_in_vector )

            if self.population.recombine_in_vector > rand and changed:
                self.population.recombine(vector)


        return changed





    def infectHost(self, host):
        """ Infects given host. """

        changed = False
        for genome,fitness in self.pathogens.items():
            if genome not in host.pathogens.keys() and genome not in host.protection_sequences and np.random.binomial( self.population.inoculum_host, self.population.inoculation_rate_host * fitness / self.sum_fitness, 1 ) > 0:
                host.pathogens[genome] = host.population.fitnessHost(genome)
                host.sum_fitness += host.pathogens[genome]
                host.population.num_infected_hosts += 1
                changed = True

            rand = np.random.random()
            if self.population.mutate_in_host > rand and changed:
                self.population.mutate(host)
                rand = rand / self.population.mutate_in_host
            else:
                rand = ( rand - self.population.mutate_in_host ) / ( 1 - self.population.mutate_in_host )

            if self.population.recombine_in_host > rand and changed:
                self.population.recombine(host)


        return changed


    def recover(self):
        ''' Remove all infections '''
        self.pathogens = {}
        self.population.num_infected_hosts -= 1


    def applyTreatment(self, treatment_seqs):
        """ Applies given treatment on this host """

        genomes_removed = []
        for genome in self.pathogens:
            for seq in treatment_seqs:
                if seq not in genome:
                    genomes_removed.append(genome)
                    break



        for g in genomes_removed:
            del self.pathogens[g]

        if len(self.pathogens) == 0:
            self.population.num_infected_hosts -= 1




class Vector(object):
    """docstring for Vector."""

    def __init__(self, population, id):
        super(Vector, self).__init__()
        self.id = id
        self.pathogens = {}
        self.population = population
        self.sum_fitness = 0
        self.protection_sequences = []


    def infectHost(self, host):
        """ Infects given host """

        changed = False
        for genome,fitness in self.pathogens.items():
            if genome not in host.pathogens.keys() and genome not in host.protection_sequences and np.random.binomial( self.population.inoculum_vector, self.population.inoculation_rate_vector * fitness / self.sum_fitness, 1 ) > 0:
                host.pathogens[genome] = host.population.fitnessHost(genome)
                host.sum_fitness += host.pathogens[genome]
                host.population.num_infected_hosts += 1
                changed = True

            rand = np.random.random()
            if self.population.mutate_in_host > rand and changed:
                self.population.mutate(host)
                rand = rand / self.population.mutate_in_host
            else:
                rand = ( rand - self.population.mutate_in_host ) / ( 1 - self.population.mutate_in_host )

            if self.population.recombine_in_host > rand and changed:
                self.population.recombine(host)


        return changed



    def recover(self):
        ''' Remove all infections '''
        self.pathogens = {}
        self.population.num_infected_vectors -= 1

    def applyTreatment(self, treatment_seqs):
        """ Applies given treatment on this vector """

        genomes_removed = []
        for genome in self.pathogens:
            for seq in treatment_seqs:
                if seq not in genome:
                    genomes_removed.append(genome)
                    break



        for g in genomes_removed:
            del self.pathogens[g]

        if len(self.pathogens) == 0:
            self.population.num_infected_vectors -= 1




class Population(object):
    """docstring for Population."""

    def __init__(self, model, params, id):

        super(Population, self).__init__()

        self.id = id

        self.hosts = [ Host(self,id) for id in range(params.num_hosts) ]
        self.vectors = [ Vector(self,id) for id in range(params.num_vectors) ]
        self.neighbors = {}
        self.total_migration_rate = 0
        self.num_infected_hosts = 0
        self.num_infected_vectors = 0

        self.num_loci = params.num_loci
        self.possible_alleles = params.possible_alleles

        self.num_hosts = params.num_hosts
        self.num_vectors = params.num_vectors
        self.fitnessHost = params.fitnessHost
        self.fitnessVector = params.fitnessVector
        self.contact_rate_host_vector = params.contact_rate_host_vector
        self.contact_rate_host_host = params.contact_rate_host_host
            # contact rate assumes fixed area--large populations are dense
            # populations, so contact scales linearly with both host and vector
            # populations. If you don't want this to happen, modify the population's
            # contact rate accordingly.
        self.inoculum_host = params.inoculum_host
        self.inoculum_vector = params.inoculum_vector
        self.inoculation_rate_host = params.inoculation_rate_host
        self.inoculation_rate_vector = params.inoculation_rate_vector
        self.recovery_rate_host = params.recovery_rate_host
        self.recovery_rate_vector = params.recovery_rate_vector
        self.recombine_in_host = params.recombine_in_host
        self.recombine_in_vector = params.recombine_in_vector
        self.mutate_in_host = params.mutate_in_host
        self.mutate_in_vector = params.mutate_in_vector
        self.vector_borne = params.vector_borne
        self.host_host_transmission = params.host_host_transmission


    def addHosts(self, num_hosts):
        """ Add a number of healthy hosts to population """

        self.hosts += [ Host( self, len(self.hosts) + i ) for i in range(num_hosts) ]


    def addVectors(self, num_vectors):
        """ Add a number of healthy vectors to population """

        self.vectors += [ Vector( self, len(self.vectors) + i ) for i in range(num_vectors) ]


    def removeHosts(self, num_hosts):
        """ Remove a number of random hosts from population """
        for _ in range(num_hosts):
            host_removed = np.random.choice(self.hosts)
            if len(host_removed.pathogens):
                self.num_infected_hosts -= 1

            self.hosts.remove( host_removed )


    def removeVectors(self, num_vectors):
        """ Remove a number of random vectors from population """
        for _ in range(num_vectors):
            vector_removed = np.random.choice(self.vectors)
            if len(vector_removed.pathogens):
                self.num_infected_vectors -= 1

            self.vectors.remove( vector_removed )


    def addPathogens(self, genomes_numbers, hosts=True):
        """ Seeds pathogens according to strains dict (keys=genomes,
            values=num of infections); seeds on hosts unless hosts=False """

        for genome in genomes_numbers:
            if len(genome) == self.num_loci and all( [ allele in self.possible_alleles for allele in genome ] ):
                if hosts:
                    self.num_infected_hosts += 1
                    for _ in range( genomes_numbers[genome] ):
                        new_fitness = self.fitnessHost(genome)
                        rand_host = np.random.choice(self.hosts)
                        rand_host.pathogens[genome] = new_fitness
                        rand_host.sum_fitness += new_fitness


                else:
                    self.num_infected_vectors += 1
                    for _ in range( genomes_numbers[genome] ):
                        new_fitness = self.fitnessVector(genome)
                        rand_vector = np.random.choice(self.vectors)
                        rand_vector.pathogens[genome] = new_fitness
                        rand_vector.sum_fitness += new_fitness



            else:
                raise ValueError('Genome ' + genome + ' must be of length ' + str(self.num_loci) + ' and contain only ' + self.possible_alleles + ' characters.')




    def recoverHost(self, index_host):
        """ Treats host with this specific index. """

        self.hosts[index_host].recover()


    def recoverVector(self, index_vector):
        """ Treats vector with this specific index. """

        self.vectors[index_vector].recover()


    def treatHosts(self, frac_hosts, treatment_seqs):
        """ Treat random hosts """

        infected_hosts = []
        for host in self.hosts:
            if len( host.pathogens ):
                infected_hosts.append( host )


        treat_hosts = np.random.choice( infected_hosts, int( frac_hosts * len( infected_hosts ) ) )
        for host in treat_hosts:
            self.hosts[host].applyTreatment(treatment_seqs)



    def treatVectors(self, frac_vectors, treatment_seqs):
        """ Treat random vectors """

        infected_vectors = []
        for vector in self.vectors:
            if len( vector.pathogens ):
                infected_vectors.append( vector )


        treat_vectors = np.random.choice( infected_vectors, int( frac_vectors * len( infected_vectors ) ) )
        for vector in treat_vectors:
            self.vectors[vector].applyTreatment(treatment_seqs)



    def protectHosts(self, frac_hosts, protection_sequence):
        """ Treat random hosts """

        protect_hosts = np.random.choice( self.hosts, int( frac_hosts * len( self.hosts ) ) )
        for host in protect_hosts:
            self.hosts[host].protection_sequences.append(protection_sequence)



    def protectVectors(self, frac_vectors, protection_sequence):
        """ Treat random vectors """

        protect_vectors = np.random.choice( self.vectors, int( frac_vectors * len( self.vectors ) ) )
        for vector in protect_vectors:
            self.vectors[vector].protection_sequences.append(protection_sequence)



    def setNeighbor(self, neighbor, rate):
        """ Adds or edits a neighbor to this population and associates the
            corresponding migration rate (from this population to the neighboring one). """

        if neighbor in self.neighbors:
            self.total_migration_rate -= self.neighbors[neighbor]

        self.neighbors[neighbor] = rate
        self.total_migration_rate += rate

    def migrate(self, target_pop, num_hosts, num_vectors):
        """ Transfers hosts and/or vectors to a target population """

        migrating_hosts = np.random.choice(self.hosts,num_hosts)
        for host in migrating_hosts:
            self.hosts.remove(host)
            target_pop.hosts.append(host)
            if len(host.pathogens):
                self.num_infected_hosts -= 1
                target_pop.num_infected_hosts += 1


        migrating_vectors = np.random.choice(self.vectors,num_vectors)
        for vector in migrating_vectors:
            self.vectors.remove(vector)
            target_pop.vectors.append(vector)
            if len(vector.pathogens):
                self.num_infected_vectors -= 1
                target_pop.num_infected_vectors += 1




    def contactVectorBorne(self, index_host, index_vector):
        """ Carries out a contact and possible transmission event between this host
            and vector. """

        temp_host = cp.deepcopy(self.hosts[index_host])
        self.vectors[index_vector].infectHost(self.hosts[index_host])
        temp_host.infectVector(self.vectors[index_vector])


    def contactHostHost(self, index_host1, index_host2):
        """ Carries out a contact and possible transmission event between two
            hosts if not vector-borne. """

        temp_host = cp.deepcopy(self.hosts[index_host1])
        self.hosts[index_host2].infectHost(self.hosts[index_host1])
        temp_host.infectHost(self.hosts[index_host2])


    def mutate(self, host_or_pathogen):
        """ Creates a new genotype from a de novo mutation event in the host or
            pathogen given. """

        if len(host_or_pathogen.pathogens) > 0:
            old_genome = np.random.choice( list( host_or_pathogen.pathogens.keys() ) )
            mut_index = np.random.randint( self.num_loci )
            new_genome = old_genome[0:mut_index] + np.random.choice( list(self.possible_alleles) ) + old_genome[mut_index+1:]
            host_or_pathogen.pathogens[new_genome] = self.fitnessHost(new_genome)



    def recombine(self, host_or_pathogen):
        """ Creates all new genotypes from all possible recombination events in
            the host or pathogen given. """

        if len(host_or_pathogen.pathogens) > 0:
            new_genomes = [""]
            for position in range( self.num_loci ):
                new_genomes_position = new_genomes
                new_genomes = []
                alleles_at_locus = []
                for genome in host_or_pathogen.pathogens:
                    if genome[position] not in alleles_at_locus:
                        alleles_at_locus.append(genome[position])
                        for new_genome in new_genomes_position:
                            new_genomes.append( new_genome + genome[position] )




            for new_genome in new_genomes:
                if new_genome not in host_or_pathogen.pathogens.keys():
                    host_or_pathogen.pathogens[new_genome] = self.fitnessVector(new_genome)









class Intervention(object):
    """docstring for Intervention."""

    def __init__(self, time, function, args):
        super(Intervention, self).__init__()
        self.time = time
        self.intervention = function
        self.args = args

    def doIntervention(self):
        """ Intervention. """

        self.intervention(*args)



class Parameters(object):
    """docstring for Parameters."""

    def __init__(self,
        num_loci, possible_alleles,
        num_hosts, num_vectors, fitnessHost, fitnessVector,
        contact_rate_host_vector, contact_rate_host_host,
        inoculum_host, inoculum_vector, inoculation_rate_host, inoculation_rate_vector,
        recovery_rate_host, recovery_rate_vector,
        recombine_in_host, recombine_in_vector,
        mutate_in_host, mutate_in_vector,
        vector_borne, host_host_transmission):

        super(Parameters, self).__init__()
        self.num_loci = num_loci
        self.possible_alleles = possible_alleles

        self.num_hosts = num_hosts
        self.num_vectors = num_vectors
        self.fitnessHost = fitnessHost
        self.fitnessVector = fitnessVector
        self.contact_rate_host_vector = contact_rate_host_vector
        self.contact_rate_host_host = contact_rate_host_host
            # contact rate assumes fixed area--large populations are dense
            # populations, so contact scales linearly with both host and vector
            # populations. If you don't want this to happen, modify the population's
            # contact rate accordingly.
        self.inoculum_host = inoculum_host
        self.inoculum_vector = inoculum_vector
        self.inoculation_rate_host = inoculation_rate_host
        self.inoculation_rate_vector = inoculation_rate_vector
        self.recovery_rate_host = recovery_rate_host
        self.recovery_rate_vector = recovery_rate_vector

        self.recombine_in_host = recombine_in_host
        self.recombine_in_vector = recombine_in_vector
        self.mutate_in_host = mutate_in_host
        self.mutate_in_vector = mutate_in_vector

        self.vector_borne = vector_borne
        self.host_host_transmission = host_host_transmission
