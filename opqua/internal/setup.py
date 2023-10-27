
"""Contains class Intervention."""

class Setup(object):
    """Class defines a setup with population parameters.
    
    Attributes:
        id (String): key of the Setup inside model dictionary.
        num_loci (int>0): length of each pathogen genome string.
        possible_alleles (String or list of Strings with num_loci elements): set of possible 
            characters in all genome string, or at each position in genome string. 
        fitnessHost (callable, takes a String argument and returns a number >= 0): function 
            that evaluates relative fitness in head-to-head competition for different genomes 
            within the same host.
        contactHost (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying probability of a given host being chosen to be the 
            infector in a contact event, based on genome sequence of pathogen.
        receiveContactHost (callable, takes a String argument and returns a number 0-1): function 
            that returns coefficient modifying probability of a given host being chosen to be 
            the infected in a contact event, based on genome sequence of pathogen.
        mortalityHost (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying death rate for a given host, based on genome sequence 
            of pathogen.
        natalityHost (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying birth rate for a given host, based on genome sequence 
            of pathogen.
        recoveryHost (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying recovery rate for a given host based on genome sequence 
            of pathogen.
        migrationHost (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying migration rate for a given host based on genome sequence 
            of pathogen.
        populationContactHost (callable, takes a String argument and returns a number 0-1): function 
            that returns coefficient modifying population contact rate for a given host based on 
            genome sequence of pathogen.
        mutationHost (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying mutation rate for a given host based on genome sequence 
            of pathogen.
        recombinationHost (callable, takes a String argument and returns a number 0-1): function 
            that returns coefficient modifying recombination rate for a given host based on genome 
            sequence of pathogen.
        fitnessVector (callable, takes a String argument and returns a number >=0): function that 
            evaluates relative fitness in head-to-head competition for different genomes within 
            the same vector.
        contactVector (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying probability of a given vector being chosen to be the 
            infector in a contact event, based on genome sequence of pathogen.
        receiveContactVector (callable, takes a String argument and returns a number 0-1): function 
            that returns coefficient modifying probability of a given vector being chosen to be the 
            infected in a contact event, based on genome sequence of pathogen.
        mortalityVector (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying death rate for a given vector, based on genome sequence 
            of pathogen.
        natalityVector (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying birth rate for a given vector, based on genome sequence 
            of pathogen.
        recoveryVector (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying recovery rate for a given vector based on genome sequence 
            of pathogen.
        migrationVector (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying migration rate for a given vector based on genome sequence 
            of pathogen.
        populationContactVector (callable, takes a String argument and returns a number 0-1): function 
            that returns coefficient modifying population contact rate for a given vector based on 
            genome sequence of pathogen.
        mutationVector (callable, takes a String argument and returns a number 0-1): function that 
            returns coefficient modifying mutation rate for a given vector based on genome sequence 
            of pathogen.
        recombinationVector (callable, takes a String argument and returns a number 0-1): function 
            that returns coefficient modifying recombination rate for a given vector based on genome 
            sequence of pathogen.
        contact_rate_host_vector (number >= 0): rate of host-vector contact events, not necessarily 
            transmission, assumes constant population density; evts/time.
        transmission_efficiency_host_vector (float): fraction of host-vector contacts
            that result in successful transmission.
        transmission_efficiency_vector_host (float): fraction of vector-host contacts
            that result in successful transmission.
        contact_rate_host_host (number >= 0): rate of host-host contact events, not
            necessarily transmission, assumes constant population density; evts/time.
        transmission_efficiency_host_host (float): fraction of host-host contacts
                that result in successful transmission.
        mean_inoculum_host (int >= 0): mean number of pathogens that are transmitted from
            a vector or host into a new host during a contact event.
        mean_inoculum_vector (int >= 0) mean number of pathogens that are transmitted
            from a host to a vector during a contact event.
        recovery_rate_host (number >= 0): rate at which hosts clear all pathogens;
            1/time.
        recovery_rate_vector (number >= 0): rate at which vectors clear all pathogens
            1/time.
        recovery_rate_vector (number >= 0): rate at which vectors clear all pathogens
            1/time.
        mortality_rate_host (number 0-1): rate at which infected hosts die from disease.
        mortality_rate_vector (number 0-1): rate at which infected vectors die from
            disease.
        recombine_in_host (number >= 0): rate at which recombination occurs in host;
            evts/time.
        recombine_in_vector (number >= 0): rate at which recombination occurs in vector;
            evts/time.
        num_crossover_host (number >= 0): mean of a Poisson distribution modeling the number
            of crossover events of host recombination events.
        num_crossover_vector (number >= 0): mean of a Poisson distribution modeling the
            number of crossover events of vector recombination events.
        mutate_in_host (number >= 0): rate at which mutation occurs in host; evts/time.
        mutate_in_vector (number >= 0): rate at which mutation occurs in vector; evts/time.
        death_rate_host (number >= 0): natural host death rate; 1/time.
        death_rate_vector (number >= 0): natural vector death rate; 1/time.
        birth_rate_host (number >= 0): infected host birth rate; 1/time.
        birth_rate_vector (number >= 0): infected vector birth rate; 1/time.
        vertical_transmission_host (number 0-1): probability that a host is infected by its
            parent at birth.
        vertical_transmission_vector (number 0-1): probability that a vector is infected by
            its parent at birth.
        inherit_protection_host (number 0-1): probability that a host inherits all
            protection sequences from its parent.
        inherit_protection_vector (number 0-1): probability that a vector inherits all
            protection sequences from its parent.
        protection_upon_recovery_host (None or array-like of length 2 with int 0-num_loci): defines 
            indexes in genome string that define substring to be added to host protection sequences 
            after recovery.
        protection_upon_recovery_vector (None or array-like of length 2 with int 0-num_loci): defines 
            indexes in genome string that define substring to be added to vector protection sequences 
            after recovery.
    """

    def __init__(
            self,
            id,
            num_loci, possible_alleles,
            fitnessHost, contactHost, receiveContactHost, mortalityHost,
            natalityHost, recoveryHost, migrationHost,
            populationContactHost, receivePopulationContactHost,
            mutationHost, recombinationHost,
            fitnessVector, contactVector, receiveContactVector, mortalityVector,
            natalityVector,recoveryVector, migrationVector,
            populationContactVector, receivePopulationContactVector,
            mutationVector, recombinationVector,
            contact_rate_host_vector,
            transmission_efficiency_host_vector,
            transmission_efficiency_vector_host,
            contact_rate_host_host,
            transmission_efficiency_host_host,
            mean_inoculum_host, mean_inoculum_vector,
            recovery_rate_host, recovery_rate_vector,
            mortality_rate_host,mortality_rate_vector,
            recombine_in_host, recombine_in_vector,
            num_crossover_host, num_crossover_vector,
            mutate_in_host, mutate_in_vector, death_rate_host,death_rate_vector,
            birth_rate_host, birth_rate_vector,
            vertical_transmission_host, vertical_transmission_vector,
            inherit_protection_host, inherit_protection_vector,
            protection_upon_recovery_host, protection_upon_recovery_vector):
        """Create a new Setup.

        Arguments:
            id (String): key of the Setup inside model dictionary.
            num_loci (int>0): length of each pathogen genome string.
            possible_alleles (String or list of Strings with num_loci elements): set of possible 
                characters in all genome string, or at each position in genome string. 
            fitnessHost (callable, takes a String argument and returns a number >= 0): function 
                that evaluates relative fitness in head-to-head competition for different genomes 
                within the same host.
            contactHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying probability of a given host being chosen to be the 
                infector in a contact event, based on genome sequence of pathogen.
            receiveContactHost (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying probability of a given host being chosen to be 
                the infected in a contact event, based on genome sequence of pathogen.
            mortalityHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying death rate for a given host, based on genome sequence 
                of pathogen.
            natalityHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying birth rate for a given host, based on genome sequence 
                of pathogen.
            recoveryHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying recovery rate for a given host based on genome sequence 
                of pathogen.
            migrationHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying migration rate for a given host based on genome sequence 
                of pathogen.
            populationContactHost (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying population contact rate for a given host based on 
                genome sequence of pathogen.
            mutationHost (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying mutation rate for a given host based on genome sequence 
                of pathogen.
            recombinationHost (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying recombination rate for a given host based on genome 
                sequence of pathogen.
            fitnessVector (callable, takes a String argument and returns a number >=0): function that 
                evaluates relative fitness in head-to-head competition for different genomes within 
                the same vector.
            contactVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying probability of a given vector being chosen to be the 
                infector in a contact event, based on genome sequence of pathogen.
            receiveContactVector (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying probability of a given vector being chosen to be the 
                infected in a contact event, based on genome sequence of pathogen.
            mortalityVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying death rate for a given vector, based on genome sequence 
                of pathogen.
            natalityVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying birth rate for a given vector, based on genome sequence 
                of pathogen.
            recoveryVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying recovery rate for a given vector based on genome sequence 
                of pathogen.
            migrationVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying migration rate for a given vector based on genome sequence 
                of pathogen.
            populationContactVector (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying population contact rate for a given vector based on 
                genome sequence of pathogen.
            mutationVector (callable, takes a String argument and returns a number 0-1): function that 
                returns coefficient modifying mutation rate for a given vector based on genome sequence 
                of pathogen.
            recombinationVector (callable, takes a String argument and returns a number 0-1): function 
                that returns coefficient modifying recombination rate for a given vector based on genome 
                sequence of pathogen.
            contact_rate_host_vector (number >= 0): rate of host-vector contact events, not necessarily 
                transmission, assumes constant population density; evts/time.
            transmission_efficiency_host_vector (float): fraction of host-vector contacts
                that result in successful transmission.
            transmission_efficiency_vector_host (float): fraction of vector-host contacts
                that result in successful transmission.
            contact_rate_host_host (number >= 0): rate of host-host contact events, not
                necessarily transmission, assumes constant population density; evts/time.
            transmission_efficiency_host_host (float): fraction of host-host contacts
                    that result in successful transmission.
            mean_inoculum_host (int >= 0): mean number of pathogens that are transmitted from
                a vector or host into a new host during a contact event.
            mean_inoculum_vector (int >= 0) mean number of pathogens that are transmitted
                from a host to a vector during a contact event.
            recovery_rate_host (number >= 0): rate at which hosts clear all pathogens;
                1/time.
            recovery_rate_vector (number >= 0): rate at which vectors clear all pathogens
                1/time.
            recovery_rate_vector (number >= 0): rate at which vectors clear all pathogens
                1/time.
            mortality_rate_host (number 0-1): rate at which infected hosts die from disease.
            mortality_rate_vector (number 0-1): rate at which infected vectors die from
                disease.
            recombine_in_host (number >= 0): rate at which recombination occurs in host;
                evts/time.
            recombine_in_vector (number >= 0): rate at which recombination occurs in vector;
                evts/time.
            num_crossover_host (number >= 0): mean of a Poisson distribution modeling the number
                of crossover events of host recombination events.
            num_crossover_vector (number >= 0): mean of a Poisson distribution modeling the
                number of crossover events of vector recombination events.
            mutate_in_host (number >= 0): rate at which mutation occurs in host; evts/time.
            mutate_in_vector (number >= 0): rate at which mutation occurs in vector; evts/time.
            death_rate_host (number >= 0): natural host death rate; 1/time.
            death_rate_vector (number >= 0): natural vector death rate; 1/time.
            birth_rate_host (number >= 0): infected host birth rate; 1/time.
            birth_rate_vector (number >= 0): infected vector birth rate; 1/time.
            vertical_transmission_host (number 0-1): probability that a host is infected by its
                parent at birth.
            vertical_transmission_vector (number 0-1): probability that a vector is infected by
                its parent at birth.
            inherit_protection_host (number 0-1): probability that a host inherits all
                protection sequences from its parent.
            inherit_protection_vector (number 0-1): probability that a vector inherits all
                protection sequences from its parent.
            protection_upon_recovery_host (None or array-like of length 2 with int 0-num_loci): defines 
                indexes in genome string that define substring to be added to host protection sequences 
                after recovery.
            protection_upon_recovery_vector (None or array-like of length 2 with int 0-num_loci): defines 
                indexes in genome string that define substring to be added to vector protection sequences 
                after recovery.
        """

        super(Setup, self).__init__()

        self.id = id

        self.num_loci = num_loci
        if isinstance(possible_alleles, list):
            self.possible_alleles = possible_alleles
        else:
            self.possible_alleles = [possible_alleles] * self.num_loci
                # possible_alleles must be a list with all available alleles for
                # each position

        self.fitnessHost = fitnessHost
        self.contactHost = contactHost
        self.receiveContactHost = receiveContactHost
        self.mortalityHost = mortalityHost
        self.natalityHost = natalityHost
        self.recoveryHost = recoveryHost
        self.migrationHost = migrationHost
        self.populationContactHost = populationContactHost
        self.receivePopulationContactHost = receivePopulationContactHost
        self.mutationHost = mutationHost
        self.recombinationHost = recombinationHost

        self.fitnessVector = fitnessVector
        self.contactVector = contactVector
        self.receiveContactVector = receiveContactVector
        self.mortalityVector = mortalityVector
        self.natalityVector = natalityVector
        self.recoveryVector = recoveryVector
        self.migrationVector = migrationVector
        self.populationContactVector = populationContactVector
        self.receivePopulationContactVector = receivePopulationContactVector
        self.mutationVector = mutationVector
        self.recombinationVector = recombinationVector

        self.contact_rate_host_vector = contact_rate_host_vector
        self.contact_rate_host_host = contact_rate_host_host
            # contact rates assumes scaling area--large populations are equally
            # dense as small ones, so contact is constant with both host and
            # vector populations. If you don't want this to happen, modify the
            # population's contact rate accordingly.
            # Examines contacts between infected hosts and all hosts
        self.transmission_efficiency_host_vector = transmission_efficiency_host_vector
        self.transmission_efficiency_vector_host = transmission_efficiency_vector_host
        self.transmission_efficiency_host_host = transmission_efficiency_host_host
        self.mean_inoculum_host = mean_inoculum_host
        self.mean_inoculum_vector = mean_inoculum_vector
        self.recovery_rate_host = recovery_rate_host
        self.recovery_rate_vector = recovery_rate_vector
        self.mortality_rate_host = mortality_rate_host
        self.mortality_rate_vector = mortality_rate_vector

        self.recombine_in_host = recombine_in_host
        self.recombine_in_vector = recombine_in_vector
        self.num_crossover_host = num_crossover_host
        self.num_crossover_vector = num_crossover_vector
        self.mutate_in_host = mutate_in_host
        self.mutate_in_vector = mutate_in_vector

        self.death_rate_host = death_rate_host
        self.death_rate_vector = death_rate_vector
        self.birth_rate_host = birth_rate_host
        self.birth_rate_vector = birth_rate_vector

        self.vertical_transmission_host = vertical_transmission_host
        self.vertical_transmission_vector = vertical_transmission_vector
        self.inherit_protection_host = inherit_protection_host
        self.inherit_protection_vector = inherit_protection_vector

        self.protection_upon_recovery_host = protection_upon_recovery_host
        self.protection_upon_recovery_vector = protection_upon_recovery_vector
