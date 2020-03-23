
import numpy as np

class Treatment(object):
    """docstring for Treatment."""

    def __init__(self, penalty_sus, penalty_res, resistance):
        super(Treatment, self).__init__()
        self.penalty = [penalty_sus, penalty_res]
        self.resistance = resistance


class Host(object):
    """docstring for Host."""

    def __init__(self, population):
        super(Host, self).__init__()
        self.pathogens = {}
        self.population = population
        self.sum_fitness = 0


    def infectVector(vector):
        """ Infects given vector """

        for genome,fitness in self.pathogens:
            if genome not in vector.pathogens.keys() and np.random.binomial( self.inoculum_host, self.inoculation_rate_host * fitness / self.sum_fitness, 1 ) > 0:
                vector.pathogens[genome] = vector.population.fitnessVector(genome)
                vector.sum_fitness += vector.pathogens[genome]





    def recover():
        ''' Remove all infections '''
        self.infections = {}


    def mutate():
    """ Creates a new genotype from a de novo mutation event. """
        old_genome = np.random.choice( self.infections.keys() )
        mut_index = np.random.randint( population.genome_length )
        new_genome = old_genome[0:mut_index] + np.random.choice( list(population.genome_alphabet) ) + old_genome[mut_index+1:-1]
        self.pathogens[new_genome] = self.population.fitnessHost(new_genome)


    def applyTreatment(treatment):
        """ Applies given treatment on this host """

        for genome,fitness in self.pathogens:
            resistant = True
            for seq in treatment.resistance:
                if seq not in genome:
                    resistant = False
                    break


            fitness = fitness * treatment.penalty[resistant]




class Vector(object):
    """docstring for Vector."""

    def __init__(self, population):
        super(Vector, self).__init__()
        self.pathogens = {}
        self.population = population
        self.sum_fitness = 0


    def infectHost(host):
        """ Infects given host """

        for genome,fitness in self.pathogens:
            if genome not in vector.pathogens.keys() and np.random.binomial( self.inoculum_vector, population.inoculation_rate_vector * fitness / self.sum_fitness, 1 ) > 0:
                host.pathogens[genome] = host.population.fitnessHost(genome)
                host.sum_fitness += host.pathogens[genome]



    def recover():
        ''' Remove all infections '''
        self.infections = {}


    def recombine():
        """ Creates all new genotypes from all possible recombination events. """

        new_genomes = [""]
        for position in range( population.genome_length ):
            new_genomes_position = new_genomes
            new_genomes = []
            alleles_at_locus = []
            for genome in self.pathogens:
                if genome[position] not in alleles_at_locus:
                    alleles_at_locus.append(genome[position])
                    for new_genome in new_genomes_position:
                        new_genomes.append( new_genome + genome[position] )




        for new_genome in new_genomes:
            if new_genome not in self.pathogens.keys():
                self.pathogens[new_genome] = self.population.fitnessVector(new_genome)






class Population(object):
    """docstring for Population."""

    def __init__(self, model, params):

        super(Population, self).__init__()
        self.hosts = [ Host(self) for _ in range(params.num_hosts) ]
        self.vectors = [ Vector(self) for _ in range(params.num_vectors) ]
        self.neighbors = {}
        self.total_migration_rate = 0

        self.genome_length = params.genome_length
        self.genome_alphabet = params.genome_alphabet
        self.treatments = params.treatments
        self.treatment_rate = params.treatment_rate
        self.contact_rate = params.contact_rate

        self.fitnessHost = params.fitnessHost
        self.fitnessVector = params.fitnessVector
        self.inoculum_host = params.inoculum_host
        self.inoculation_rate_host = params.inoculation_rate_host
        self.inoculation_rate_vector = params.inoculation_rate_vector
        self.recovery_rate_host = params.recovery_rate_host
        self.recovery_rate_vector = params.recovery_rate_vector


    def addHosts(num_hosts):
        """ Add a number of healthy hosts to population """

        self.hosts += [ Host(self) for _ in range(num_hosts) ]


    def addVectors(num_vectors):
        """ Add a number of healthy vectors to population """

        self.vectors += [ Vector(self) for _ in range(num_vectors) ]

    def removeHosts(num_hosts):
        """ Remove a number of random hosts from population """
        for _ in range(num_hosts):
            self.hosts.remove( np.random.choice(self.hosts) )


    def removeVectors(num_vectors):
        """ Remove a number of random vectors from population """
        for _ in range(num_vectors):
            self.vectors.remove( np.random.choice(self.vectors) )


    def addpathogens( strains, hosts=True ):
        """ Seeds pathogens according to strains dict (keys=genomes,
            values=num of infections); seeds on hosts unless hosts=False """

        for genome in strains:
            if hosts:
                for _ in range( strains[genome] ):
                    np.random.choice(self.hosts).pathogens[genome] = self.fitnessHost(genome)


            else:
                for _ in range( strains[genome] ):
                    np.random.choice(self.vectors).pathogens[genome] = self.fitnessVector(genome)




    def treatHosts(frac_hosts):
        """ Treat random hosts """

        infected_hosts = []
        for host in self.hosts:
            if len( host.pathogens ):
                infected_hosts.append( host )


        treat_hosts = np.random.choice( infected_hosts, int( frac_hosts * len( infected_hosts ) ) )
        for host in migrating_hosts:
            self.hosts.remove(host)



    def setNeighbor(neighbor,rate):
        """ Adds or edits a neighbor to this population and associates the
            corresponding migration rate (from this population to the neighboring one). """

        if neighbor in self.neighbors.keys():
            self.total_migration_rate -= self.neighbors[neighbor]

        self.neighbors[neighbor] = rate
        self.total_migration_rate += rate


    def migrate(target_pop,num_hosts,num_vectors):
        """ Transfers hosts and/or vectors to a target population """

        migrating_hosts = np.random.choice(self.hosts,num_hosts)
        migrating_vectors = np.random.choice(self.vectors,num_vectors)
        for host in migrating_hosts:
            self.hosts.remove(host)
            target_pop.hosts.append(host)

        for vector in migrating_vectors:
            self.vectors.remove(vector)
            target_pop.vectors.append(vector)





class Parameters(object):
    """docstring for Parameters."""

    def __init__(self,
        genome_length, genome_alphabet, treatments, contact_rate, treatment_rate,
        num_hosts, num_vectors, fitnessHost, fitnessVector,
        inoculum_host, inoculum_vector, inoculation_rate_host, inoculation_rate_vector,
        recovery_rate_host, recovery_rate_vector):

        super(Parameters, self).__init__()
        self.genome_length = genome_length
        self.genome_alphabet = genome_alphabet
        self.treatments = treatments
        self.treatment_rate = treatment_rate
        self.contact_rate = contact_rate

        self.fitnessHost = fitnessHost
        self.fitnessVector = fitnessVector
        self.inoculum_host = inoculum_host
        self.inoculation_rate_host = inoculation_rate_host
        self.inoculation_rate_vector = inoculation_rate_vector
        self.recovery_rate_host = recovery_rate_host
        self.recovery_rate_vector = recovery_rate_vector
