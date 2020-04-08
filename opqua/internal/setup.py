
"""Contains class Intervention."""

class Setup(object):
    """Class defines a setup with population parameters."""

    def __init__(
            self,
            num_loci, possible_alleles,
            fitnessHost, fitnessVector,
            contact_rate_host_vector, contact_rate_host_host,
            mean_inoculum_host, mean_inoculum_vector,
            recovery_rate_host, recovery_rate_vector,
            recombine_in_host, recombine_in_vector,
            mutate_in_host, mutate_in_vector, death_rate_host,death_rate_vector,
            protection_upon_recovery_host, protection_upon_recovery_vector):
        """Create a new Setup.

        Arguments:
        num_loci -- length of each pathogen genome string (int > 0)
        possible_alleles -- set of possible characters in all genome string, or
            at each position in genome string (String or list of Strings with
            num_loci elements)
        fitnessHost -- relative fitness in head-to-head competition within host
            (number >= 0)
        fitnessVector -- relative fitness in head-to-head competition within
            vector (number >= 0)
        contact_rate_host_vector -- rate of host-vector contact events, not
            necessarily transmission, assumes constant population density;
            evts/time (number >= 0)
        contact_rate_host_host -- rate of host-host contact events, not
            necessarily transmission, assumes constant population density;
            evts/time (number >= 0)
        mean_inoculum_host -- mean number of pathogens that are transmitted from
            a vector or host into a new host during a contact event (int >= 0)
        mean_inoculum_vector -- mean number of pathogens that are transmitted
            from a host to a vector during a contact event (int >= 0)
        recovery_rate_host -- rate at which hosts clear all pathogens;
            1/time (number >= 0)
        recovery_rate_vector -- rate at which vectors clear all pathogens
            1/time (number >= 0)
        recombine_in_host -- rate at which recombination occurs in host;
            evts/time (number >= 0)
        recombine_in_vector -- rate at which recombination occurs in vector;
            evts/time (number >= 0)
        mutate_in_host -- rate at which mutation occurs in host; evts/time
            (number >= 0)
        mutate_in_vector -- rate at which mutation occurs in vector; evts/time
            (number >= 0)
        death_rate_host -- infected host death rate; 1/time (number >= 0)
        death_rate_vector -- infected vector death rate; 1/time (number >= 0)
        protection_upon_recovery_host -- defines indexes in genome string that
            define substring to be added to host protection sequences after
            recovery (None or array-like of length 2 with int 0-num_loci)
        protection_upon_recovery_vector -- defines indexes in genome string that
            define substring to be added to vector protection sequences after
            recovery (None or array-like of length 2 with int 0-num_loci)
        """

        super(Setup, self).__init__()
        self.num_loci = num_loci
        if isinstance(possible_alleles, list):
            self.possible_alleles = possible_alleles
        else:
            self.possible_alleles = [possible_alleles] * self.num_loci
                # possible_alleles must be a list with all available alleles for
                # each position

        self.fitnessHost = fitnessHost
        self.fitnessVector = fitnessVector
        self.contact_rate_host_vector = contact_rate_host_vector
        self.contact_rate_host_host = contact_rate_host_host
            # contact rates assumes scaling area--large populations are equally
            # dense as small ones, so contact is constant with both host and
            # vector populations. If you don't want this to happen, modify the
            # population's contact rate accordingly.
            # Examines contacts between infected hosts and all hosts
        self.mean_inoculum_host = mean_inoculum_host
        self.mean_inoculum_vector = mean_inoculum_vector
        self.recovery_rate_host = recovery_rate_host
        self.recovery_rate_vector = recovery_rate_vector

        self.recombine_in_host = recombine_in_host
        self.recombine_in_vector = recombine_in_vector
        self.mutate_in_host = mutate_in_host
        self.mutate_in_vector = mutate_in_vector

        self.death_rate_host = death_rate_host
        self.death_rate_vector = death_rate_vector
        self.protection_upon_recovery_host = protection_upon_recovery_host
        self.protection_upon_recovery_vector = protection_upon_recovery_vector
