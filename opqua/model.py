# TODO: Update all fitnesses
# TODO: remove protections
# TODO: parallelizeable simulations
# TODO: Pathogen genome influences transmission probability, death rate
# TODO: contact between populations (without migration)
# TODO: birth/death rates in populations
# TODO: arbitrary comparments
# TODO: independent recombination of alleles
"""Contains class Model; main class user interacts with."""

import numpy as np
import pandas as pd
import textdistance as td
import seaborn as sns

from opqua.internal.host import Host
from opqua.internal.vector import Vector
from opqua.internal.population import Population
from opqua.internal.setup import Setup
from opqua.internal.intervention import Intervention
from opqua.internal.gillespie import Gillespie
from opqua.internal.data import saveToDf, getPathogens, getProtections
from opqua.internal.plot import populationsPlot, compartmentPlot, \
    compositionPlot, clustermap

class Model(object):
    """Class defines a Model.

    This is the main class that the user interacts with.

    The Model class contains populations, setups, and interventions to be used
    in simulation. Also contains groups of hosts/vectors for manipulations and
    stores model history as snapshots for each time point.

    *** --- ATTRIBUTES: --- ***
    populations -- dictionary with keys=population IDs, values=Population
        objects
    setups -- dictionary with keys=setup IDs, values=Setup objects
    interventions -- contains model interventions in the order they will occur
    groups -- dictionary with keys=group IDs, values=lists of hosts/vectors
    self.history -- dictionary with keys=time values, values=Model objects that
        are snapshots of Model at that timepoint

    *** --- METHODS: --- ***

    --- Model Initialization and Simulation: ---

    newSetup -- creates a new Setup, save it in setups dict under given name
    newIntervention -- creates a new intervention executed during simulation
    run -- simulates model for a specified length of time

    --- Data Output and Plotting: ---

    saveToDataFrame -- saves status of model to dataframe, writes to file
    getPathogens -- creates Dataframe with counts for all pathogen genomes
    getProtections -- creates Dataframe with counts for all protection sequences
    populationsPlot -- plots aggregated totals per population across time
    compartmentPlot -- plots number of naive,inf,rec,dead hosts/vectors vs time
    compositionPlot -- plots counts for pathogen genomes or resistance vs. time

    --- Model interventions: ---

    - Make and connect populations -
    newPopulation -- create a new Population object with setup parameters
    linkPopulations -- set migration rate from one population towards another
    createInterconnectedPopulations -- create new populations, link all of them
        to each other

    - Modify population parameters -
    setSetup -- assigns a given set of parameters to this population

    - Manipulate hosts and vectors in population -
    newHostGroup -- returns a list of random (healthy or any) hosts
    newVectorGroup -- returns a list of random (healthy or any) vectors
    addHosts -- adds hosts to the population
    addVectors -- adds vectors to the population
    removeHosts -- removes hosts from the population
    removeVectors -- removes vectors from the population
    addPathogensToHosts -- adds pathogens with specified genomes to hosts
    addPathogensToVectors -- adds pathogens with specified genomes to vectors
    treatHosts -- removes infections susceptible to given treatment from hosts
    treatVectors -- removes infections susceptible to treatment from vectors
    protectHosts -- adds protection sequence to hosts
    protectVectors -- adds protection sequence to vectors


    --- Preset fitness functions: ---
        * these are static methods

    stabilizingSelection -- evaluates genome fitness by decreasing with distance
        from optimal sequence
    disruptiveSelection -- evaluates genome fitness by increasing with distance
        from worst sequence
    """

    ### CONSTANTS ###
    ### Color scheme constants ###
    CB_PALETTE = ["#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]
     # www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
     # http://jfly.iam.u-tokyo.ac.jp/color/

    DEF_CMAP = sns.cubehelix_palette(
        start=.5, rot=-.75, as_cmap=True, reverse=True
        )


    ### CLASS CONSTRUCTOR ###

    def __init__(self):
        """Create a new Model object."""
        super(Model, self).__init__()
        self.populations = {}
            # dictionary with keys=population IDs, values=Population objects
        self.setups = {}
            # dictionary with keys=setup IDs, values=Setup objects
        self.interventions = []
            # contains model interventions in the order they will occur
        self.groups = {}
            # dictionary with keys=group IDs, values=lists of hosts/vectors
        self.history = {}
            # dictionary with keys=time values, values=Model objects that are
            # snapshots of Model at that timepoint


    ### MODEL METHODS ###

    ### Model initialization and simulation: ###

    def newSetup(
            self, name, preset=None,
            num_loci=None, possible_alleles=None,
            fitnessHost=None, fitnessVector=None,
            contact_rate_host_vector=None, contact_rate_host_host=None,
            mean_inoculum_host=None, mean_inoculum_vector=None,
            recovery_rate_host=None, recovery_rate_vector=None,
            recombine_in_host=None, recombine_in_vector=None,
            mutate_in_host=None, mutate_in_vector=None,
            death_rate_host=None, death_rate_vector=None,
            protection_upon_recovery_host=None,
            protection_upon_recovery_vector=None):
        """Create a new Setup, save it in setups dict under given name.

        Two preset setups exist: "vector-borne" and "host-host". You may select
        one of the preset setups with the preset keyword argument and then
        modify individual parameters with additional keyword arguments, without
        having to specify all of them.

        Preset parameter setups:

        "vector-borne":
        num_loci = num_loci or 10
        possible_alleles = possible_alleles or 'ATCG'
        fitnessHost = fitnessHost or (lambda g: 1)
        fitnessVector = fitnessVector or (lambda g: 1)
        contact_rate_host_vector = contact_rate_host_vector or 1e1
        contact_rate_host_host = contact_rate_host_host or 0
        mean_inoculum_host = mean_inoculum_host or 1e2
        mean_inoculum_vector = mean_inoculum_vector or 1e2
        recovery_rate_host = recovery_rate_host or 1e-1
        recovery_rate_vector = recovery_rate_vector or 1e-2
        recombine_in_host = recombine_in_host or 0
        recombine_in_vector = recombine_in_vector or 1e-2
        mutate_in_host = mutate_in_host or 1e-6
        mutate_in_vector = mutate_in_vector or 0
        death_rate_host = death_rate_host or 0
        death_rate_vector = death_rate_vector or 0
        protection_upon_recovery_host = ( protection_upon_recovery_host
            or None )
        protection_upon_recovery_vector = ( protection_upon_recovery_vector
            or None )

        "host-host":
        num_loci = num_loci or 10
        possible_alleles = possible_alleles or 'ATCG'
        fitnessHost = fitnessHost or (lambda g: 1)
        fitnessVector = fitnessVector or (lambda g: 1)
        contact_rate_host_vector = contact_rate_host_vector or 0
        contact_rate_host_host = contact_rate_host_host or 2e1
        mean_inoculum_host = mean_inoculum_host or 1e1
        mean_inoculum_vector = mean_inoculum_vector or 0
        recovery_rate_host = recovery_rate_host or 1e-1
        recovery_rate_vector = recovery_rate_vector or 1e1
        recombine_in_host = recombine_in_host or 1e-3
        recombine_in_vector = recombine_in_vector or 0
        mutate_in_host = mutate_in_host or 1e-6
        mutate_in_vector = mutate_in_vector or 0
        death_rate_host = death_rate_host or 0
        death_rate_vector = death_rate_vector or 0
        protection_upon_recovery_host = ( protection_upon_recovery_host
            or None )
        protection_upon_recovery_vector = ( protection_upon_recovery_vector
            or None )

        Arguments:
        name -- name of setup to be used as a key in model setups dictionary

        Keyword arguments:
        preset -- preset setup to be used: "vector-borne" or "host-host", if
            None, must define all other keyword arguments (default None; None or
            String)
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

        if preset == "vector-borne":
            num_loci = 10 if num_loci is None else num_loci
            possible_alleles = \
                'ATCG' if possible_alleles is None else possible_alleles
            fitnessHost = (lambda g: 1) if fitnessHost is None else fitnessHost
            fitnessVector = \
                (lambda g: 1) if fitnessVector is None else fitnessVector
            contact_rate_host_vector = \
                1e1 if contact_rate_host_vector is None \
                else contact_rate_host_vector
            contact_rate_host_host = \
                0 if contact_rate_host_host is None else contact_rate_host_host
            mean_inoculum_host = \
                1e2 if mean_inoculum_host is None else mean_inoculum_host
            mean_inoculum_vector = \
                1e2 if mean_inoculum_vector is None else mean_inoculum_vector
            recovery_rate_host = \
                1e-1 if recovery_rate_host is None else recovery_rate_host
            recovery_rate_vector = \
                1e-2 if recovery_rate_vector is None else recovery_rate_vector
            recombine_in_host = \
                0 if recombine_in_host is None else recombine_in_host
            recombine_in_vector = \
                1e-2 if recombine_in_vector is None else recombine_in_vector
            mutate_in_host = 1e-6 if mutate_in_host is None else mutate_in_host
            mutate_in_vector = \
                0 if mutate_in_vector is None else mutate_in_vector
            death_rate_host = 0 if death_rate_host is None else death_rate_host
            death_rate_vector = \
                0 if death_rate_vector is None else death_rate_vector
            protection_upon_recovery_host = protection_upon_recovery_host
            protection_upon_recovery_vector = protection_upon_recovery_vector

        elif preset == "host-host":
            num_loci = 10 if num_loci is None else num_loci
            possible_alleles = \
                'ATCG' if possible_alleles is None else possible_alleles
            fitnessHost = (lambda g: 1) if fitnessHost is None else fitnessHost
            fitnessVector = \
                (lambda g: 1) if fitnessVector is None else fitnessVector
            contact_rate_host_vector = \
                0 if contact_rate_host_vector is None \
                else contact_rate_host_vector
            contact_rate_host_host = \
                2e1 if contact_rate_host_host is None \
                else contact_rate_host_host
            mean_inoculum_host = \
                1e1 if mean_inoculum_host is None else mean_inoculum_host
            mean_inoculum_vector = \
                0 if mean_inoculum_vector is None else mean_inoculum_vector
            recovery_rate_host = \
                1e-1 if recovery_rate_host is None else recovery_rate_host
            recovery_rate_vector = \
                1e1 if recovery_rate_vector is None else recovery_rate_vector
            recombine_in_host = \
                1e-3 if recombine_in_host is None else recombine_in_host
            recombine_in_vector = \
                0 if recombine_in_vector is None else recombine_in_vector
            mutate_in_host = \
                1e-6 if mutate_in_host is None else mutate_in_host
            mutate_in_vector = \
                0 if mutate_in_vector is None else mutate_in_vector
            death_rate_host = \
                0 if death_rate_host is None else death_rate_host
            death_rate_vector = \
                0 if death_rate_vector is None else death_rate_vector
            protection_upon_recovery_host = protection_upon_recovery_host
            protection_upon_recovery_vector = protection_upon_recovery_vector

        self.setups[name] = Setup(
            num_loci, possible_alleles,
            fitnessHost, fitnessVector,
            contact_rate_host_vector, contact_rate_host_host,
            mean_inoculum_host, mean_inoculum_vector,
            recovery_rate_host, recovery_rate_vector,
            recombine_in_host, recombine_in_vector,
            mutate_in_host, mutate_in_vector,
            death_rate_host, death_rate_vector,
            protection_upon_recovery_host, protection_upon_recovery_vector
            )

    def newIntervention(self, time, function, args):
        """Create a new intervention to be carried out at a specific time.

        Arguments:
        time -- time at which intervention will take place (number)
        function -- intervention to be carried out (method of class Model)
        args -- contains arguments for function in positinal order (array-like)
        """

        self.interventions.append( Intervention(time, function, args) )

    def run(self,t0,tf):
        """Simulate model for a specified time between two time points.

        Simulates a time series using the Gillespie algorithm.

        Saves a dictionary containing model state history, with keys=times and
        values=Model objects with model snapshot at that time point under this
        model's history attribute.

        Arguments:
        t0 -- initial time point to start simulation at (number)
        tf -- initial time point to end simulation at (number)
        """

        sim = Gillespie(self)
        self.history = sim.run(t0,tf)


    ### Output and Plots: ###

    def saveToDataFrame(self,save_to_file,n_cores=0):
        """Save status of model to dataframe, write to file location given.

        Creates a pandas Dataframe in long format with the given model history,
        with one host or vector per simulation time in each row, and columns:
            Time - simulation time of entry
            Population - ID of this host/vector's population
            Organism - host/vector
            ID - ID of host/vector
            Pathogens - all genomes present in this host/vector separated by ;
            Protection - all genomes present in this host/vector separated by ;
            Alive - whether host/vector is alive at this time, True/False

        Arguments:
        save_to_file -- file path and name to save model data under (String)

        Keyword arguments:
        n_cores -- number of cores to parallelize file export across, if 0, all
            cores available are used (default 0; int)

        Returns:
        pandas dataframe with model history as described above
        """

        data = saveToDf(self.history,save_to_file,n_cores)

        return data

    def getPathogens(self, dat, save_to_file=""):
        """Create Dataframe with counts for all pathogen genomes in data.

        Returns sorted pandas Dataframe with counts for occurrences of all pathogen
        genomes in data passed.

        Arguments:
        data -- dataframe with model history as produced by saveToDf function

        Keyword arguments:
        save_to_file -- file path and name to save model data under, no saving
            occurs if empty string (default ''; String)

        Returns:
        pandas dataframe with Series as described above
        """

        return getPathogens(dat, save_to_file=save_to_file)

    def getProtections(self, dat, save_to_file=""):
        """Create Dataframe with counts for all protection sequences in data.

        Returns sorted pandas Dataframe with counts for occurrences of all
        protection sequences in data passed.

        Arguments:
        data -- dataframe with model history as produced by saveToDf function

        Keyword arguments:
        save_to_file -- file path and name to save model data under, no saving
            occurs if empty string (default ''; String)

        Returns:
        pandas dataframe with Series as described above
        """

        return getProtections(dat, save_to_file=save_to_file)

    def populationsPlot(
            self, file_name, data, compartment='Infected',
            hosts=True, vectors=False, num_top_populations=7,
            track_specific_populations=[], save_data_to_file="",
            x_label='Time', y_label='Hosts', figsize=(8, 4), dpi=200,
            palette=CB_PALETTE, stacked=False):
        """Create plot with aggregated totals per population across time.

        Creates a line or stacked line plot with dynamics of a compartment
        across populations in the model, with one line for each population.

        A host or vector is considered part of the recovered compartment
        if it has protection sequences of any kind and is not infected.

        Arguments:
        file_name -- file path, name, and extension to save plot under (String)
        data -- dataframe with model history as produced by saveToDf function
            (DataFrame)

        Keyword arguments:
        compartment -- subset of hosts/vectors to count totals of, can be either
            'Naive','Infected','Recovered', or 'Dead' (default 'Infected';
            String)
        hosts -- whether to count hosts (default True, Boolean)
        vectors -- whether to count vectors (default False, Boolean)
        num_top_populations -- how many populations to count separately and
            include as columns, remainder will be counted under column "Other";
             if <0, includes all populations in model (default 7; int)
        track_specific_populations -- contains IDs of specific populations to
            have as a separate column if not part of the top num_top_populations
            populations (list of Strings)
        save_data_to_file -- file path and name to save model plot data under,
            no saving occurs if empty string (default ''; String)
        x_label -- X axis title (default 'Time', String)
        y_label -- Y axis title (default 'Hosts', String)
        legend_title -- legend title (default 'Population', String)
        legend_values -- labels for each trace, if empty list, uses population
            IDs (default empty list, list of Strings)
        figsize -- dimensions of figure (default (8,4), array-like of two ints)
        dpi -- figure resolution (default 200, int)
        palette -- color palette to use for traces (default CB_PALETTE, list of
            color Strings)
        stacked -- whether to draw a regular line plot or a stacked one (default
            False, Boolean)

        Returns:
        axis object for plot with model population dynamics as described above
        """

        return populationsPlot(
            file_name, data, compartment=compartment, hosts=hosts,
            vectors=vectors, num_top_populations=num_top_populations,
            track_specific_populations=track_specific_populations,
            save_data_to_file=save_data_to_file,
            x_label=x_label, y_label=y_label, figsize=figsize, dpi=dpi,
            palette=palette, stacked=stacked
            )

    def compartmentPlot(
            self, file_name, data, populations=[], hosts=True, vectors=False,
            save_data_to_file="", x_label='Time', y_label='Hosts',
            figsize=(8, 4), dpi=200, palette=CB_PALETTE, stacked=False):
        """Create plot with number of naive,inf,rec,dead hosts/vectors vs. time.

        Creates a line or stacked line plot with dynamics of all compartments
        (naive, infected, recovered, dead) across selected populations in the
        model, with one line for each compartment.

        A host or vector is considered part of the recovered compartment
        if it has protection sequences of any kind and is not infected.

        Arguments:
        file_name -- file path, name, and extension to save plot under (String)
        data -- dataframe with model history as produced by saveToDf function
            (DataFrame)

        Keyword arguments:
        populations -- IDs of populations to include in analysis; if empty, uses
            all populations in model (default empty list; list of Strings)
        hosts -- whether to count hosts (default True, Boolean)
        vectors -- whether to count vectors (default False, Boolean)
        save_data_to_file -- file path and name to save model data under, no
            saving occurs if empty string (default ''; String)
        x_label -- X axis title (default 'Time', String)
        y_label -- Y axis title (default 'Hosts', String)
        legend_title -- legend title (default 'Population', String)
        legend_values -- labels for each trace, if empty list, uses population
            IDs (default empty list, list of Strings)
        figsize -- dimensions of figure (default (8,4), array-like of two ints)
        dpi -- figure resolution (default 200, int)
        palette -- color palette to use for traces (default CB_PALETTE, list of
            color Strings)
        stacked -- whether to draw a regular line plot or a stacked one (default
            False, Boolean)

        Returns:
        axis object for plot with model compartment dynamics as described above
        """

        return compartmentPlot(
            file_name, data, populations=populations, hosts=hosts,
            vectors=vectors, save_data_to_file=save_data_to_file,
            x_label=x_label, y_label=y_label, figsize=figsize, dpi=dpi,
            palette=palette, stacked=stacked
            )

    def compositionPlot(
            self, file_name, data, populations=[],
            type_of_composition='Pathogens', hosts=True, vectors=False,
            num_top_sequences=7, track_specific_sequences=[],
            save_data_to_file="", x_label='Time', y_label='Infections',
            figsize=(8, 4), dpi=200, palette=CB_PALETTE, stacked=True):
        """Create plot with counts for pathogen genomes or resistance vs. time.

        Creates a line or stacked line plot with dynamics of the pathogen
        strains or protection sequences across selected populations in the
        model, with one line for each pathogen genome or protection sequence
        being shown.

        Of note: sum of totals for all sequences in one time point does not
        necessarily equal the number of infected hosts and/or vectors, given
        multiple infections in the same host/vector are counted separately.

        Arguments:
        data -- dataframe with model history as produced by saveToDf function

        Keyword arguments:
        populations -- IDs of populations to include in analysis; if empty, uses
            all populations in model (default empty list; list of Strings)
        type_of_composition -- field of data to count totals of, can be either
            'Pathogens' or 'Protection' (default 'Pathogens'; String)
        hosts -- whether to count hosts (default True, Boolean)
        vectors -- whether to count vectors (default False, Boolean)
        num_top_sequences -- how many sequences to count separately and include
            as columns, remainder will be counted under column "Other"; if <0,
            includes all genomes in model (default 7; int)
        track_specific_sequences -- contains specific sequences to have
            as a separate column if not part of the top num_top_sequences
            sequences (list of Strings)
        save_data_to_file -- file path and name to save model data under, no
            saving occurs if empty string (default ''; String)
        x_label -- X axis title (default 'Time', String)
        y_label -- Y axis title (default 'Hosts', String)
        legend_title -- legend title (default 'Population', String)
        legend_values -- labels for each trace, if empty list, uses population
            IDs (default empty list, list of Strings)
        figsize -- dimensions of figure (default (8,4), array-like of two ints)
        dpi -- figure resolution (default 200, int)
        palette -- color palette to use for traces (default CB_PALETTE, list of
            color Strings)
        stacked -- whether to draw a regular line plot or a stacked one (default
            False, Boolean)

        Returns:
        axis object for plot with model sequence composition dynamics as
            described
        """

        return compositionPlot(
            file_name, data, populations=populations,
            type_of_composition=type_of_composition, hosts=hosts,
            vectors=vectors, num_top_sequences=num_top_sequences,
            track_specific_sequences=track_specific_sequences,
            save_data_to_file=save_data_to_file,
            x_label=x_label, y_label=y_label, figsize=figsize, dpi=dpi,
            palette=palette, stacked=stacked
            )

    def clustermap(
            self,
            file_name, data, num_top_sequences=-1, track_specific_sequences=[],
            seq_names=[], n_cores=0, method='weighted', metric='euclidean',
            save_data_to_file="", legend_title='Distance', legend_values=[],
            figsize=(10,10), dpi=200, color_map=DEF_CMAP):
        """Create a heatmap and dendrogram for pathogen genomes in data passed.

        Arguments:
        file_name -- file path, name, and extension to save plot under (String)
        data -- dataframe with model history as produced by saveToDf function

        Keyword arguments:
        num_top_sequences -- how many sequences to include in matrix; if <0,
            includes all genomes in data passed (default -1; int)
        track_specific_sequences -- contains specific sequences to include in
            matrixif not part of the top num_top_sequences sequences (default
            empty list; list of Strings)
        seq_names -- list with names to be used for sequence labels in matrix
            must be of same length as number of sequences to be displayed; if
            empty, uses sequences themselves (default empty list; list of
            Strings)
        n_cores -- number of cores to parallelize distance compute across, if 0,
            all cores available are used (default 0; int)
        method -- clustering algorithm to use with seaborn clustermap (default
            'weighted'; String)
        metric -- distance metric to use with seaborn clustermap (default
            'euclidean'; String)
        save_data_to_file -- file path and name to save model data under, no
            saving occurs if empty string (default ''; String)
        legend_title -- legend title (default 'Distance', String)
        figsize -- dimensions of figure (default (8,4), array-like of two ints)
        dpi -- figure resolution (default 200, int)
        color_map -- color map to use for traces (default DEF_CMAP, cmap object)

        Returns:
        figure object for plot with heatmap and dendrogram as described
        """

        return clustermap(
                file_name, data, num_top_sequences=num_top_sequences,
                track_specific_sequences=track_specific_sequences,
                seq_names=seq_names, n_cores=n_cores, method=method,
                metric=metric, save_data_to_file=save_data_to_file,
                legend_title=legend_title, legend_values=legend_values,
                figsize=figsize, dpi=dpi, color_map=color_map
                )

    ### Model interventions: ###

    def newPopulation(self, id, setup_name, num_hosts=100, num_vectors=100):
        """Create a new Population object with setup parameters.

        If population ID is already in use, appends _2 to it

        Arguments:
        id -- unique identifier for this population in the model (String)
        setup_name -- setup object with parameters for this population (Setup)

        Keyword arguments:
        num_hosts -- number of hosts to initialize population with (default 100;
            int)
        num_vectors -- number of hosts to initialize population with (default
            100; int)
        """

        if id in self.populations.keys():
            id = id+'_2'

        self.populations[id] = Population(
            self, id, self.setups[setup_name], num_hosts, num_vectors
            )

    def linkPopulations(self, pop1_id, pop2_id, rate):
        """Set migration rate from one population towards another.

        Arguments:
        neighbor -- population towards which migration rate will be specified
            (Population)
        rate -- migration rate from this population to the neighbor; evts/time
            (number)
        """

        self.populations[pop1_id].setNeighbor( self.populations[pop2_id], rate )

    def createInterconnectedPopulations(
            self, num_populations, migration_rate, id_prefix, setup_name,
            num_hosts=100, num_vectors=100):
        """Create new populations, link all of them to each other.

        All populations in this cluster are linked with the same migration rate,
        starting number of hosts and vectors, and setup parameters. Their IDs
        are numbered onto prefix given as 'id_prefix_0', 'id_prefix_1',
        'id_prefix_2', etc.

        Arguments:
        num_populations -- number of populations to be created (int)
        migration_rate -- migration rate between populations; evts/time (number)
        id_prefix -- prefix for IDs to be used for this population in the model,
            (String)
        setup_name -- setup object with parameters for all populations (Setup)

        Keyword arguments:
        num_hosts -- number of hosts to initialize population with (default 100;
            int)
        num_vectors -- number of hosts to initialize population with (default
            100; int)
        """

        new_pops = [
            Population(
                self, id_prefix + str(i), self.setups[setup_name],
                num_hosts, num_vectors
                ) for i in range(num_populations)
            ]
        new_pop_ids = []
        for pop in new_pops:
            if pop.id in self.populations.keys():
                pop.id = pop.id+'_2'

            self.populations[pop.id] = pop
            new_pop_ids.append(pop.id)

        for p1_id in new_pop_ids:
            for p2_id in new_pop_ids:
                self.linkPopulations(p1_id,p2_id,migration_rate)

    def newHostGroup(self, pop_id, group_id, num_hosts, healthy=False):
        """Return a list of random (healthy or any) hosts in population.

        Arguments:
        pop_id -- ID of population to be modified (String)
        num_vectors -- number of vectors to be sampled randomly (int)

        Keyword arguments:
        healthy -- whether to sample healthy hosts only (default True; Boolean)

        Returns:
        list containing sampled hosts
        """

        self.groups[group_id] = self.populations[pop_id].newHostGroup(num_hosts)

    def newVectorGroup(self, pop_id, group_id, num_vectors, healthy=False):
        """Return a list of random (healthy or any) vectors in population.

        Arguments:
        pop_id -- ID of population to be modified (String)
        num_vectors -- number of vectors to be sampled randomly (int)

        Keyword arguments:
        healthy -- whether to sample healthy vectors only (default True;
            Boolean)

        Returns:
        list containing sampled vectors
        """

        self.groups[group_id] = self.populations[pop_id].newVectorGroup(
            num_vectors
            )

    def addHosts(self, pop_id, num_hosts):
        """Add a number of healthy hosts to population, return list with them.

        Arguments:
        pop_id -- ID of population to be modified (String)
        num_hosts -- number of hosts to be added (int)

        Returns:
        list containing new hosts
        """

        self.populations[pop_id].addHosts(num_hosts)

    def addVectors(self, pop_id, num_vectors):
        """Add a number of healthy vectors to population, return list with them.

        Arguments:
        pop_id -- ID of population to be modified (String)
        num_vectors -- number of vectors to be added (int)

        Returns:
        list containing new vectors
        """

        self.populations[pop_id].addVectors(num_vectors)

    def removeHosts(self, pop_id, num_hosts_or_list):
        """Remove a number of specified or random hosts from population.

        Arguments:
        pop_id -- ID of population to be modified (String)
        num_hosts_or_list -- number of hosts to be sampled randomly for removal
            or list of hosts to be removed, must be hosts in this population
            (int or list of Hosts)
        """

        self.populations[pop_id].removeHosts(num_hosts_or_list)

    def removeVectors(self, pop_id, num_vectors_or_list):
        """Remove a number of specified or random vectors from population.

        Arguments:
        pop_id -- ID of population to be modified (String)
        num_vectors_or_list -- number of vectors to be sampled randomly for
            removal or list of vectors to be removed, must be vectors in this
            population (int or list of Vectors)
        """

        self.populations[pop_id].removeVectors(num_vectors_or_list)


    def addPathogensToHosts(self, pop_id, genomes_numbers, group_id=""):
        """Add specified pathogens to random hosts, optionally from a list.

        Arguments:
        pop_id -- ID of population to be modified (String)
        genomes_numbers -- dictionary containing pathogen genomes to add as keys
            and number of hosts each one will be added to as values (dict with
            keys=Strings, values=int)

        Keyword arguments:
        hosts -- list of specific hosts to sample from, if empty, samples from
            whole population (default empty list; empty)
        """

        if group_id == "":
            hosts = self.populations[pop_id].hosts
        else:
            hosts = self.groups[group_id]

        self.populations[pop_id].addPathogensToHosts(genomes_numbers,hosts)

    def addPathogensToVectors(self, pop_id, genomes_numbers, group_id=""):
        """Add specified pathogens to random vectors, optionally from a list.

        Arguments:
        pop_id -- ID of population to be modified (String)
        genomes_numbers -- dictionary containing pathogen genomes to add as keys
            and number of vectors each one will be added to as values (dict with
            keys=Strings, values=int)

        Keyword arguments:
        vectors -- list of specific vectors to sample from, if empty, samples
            from whole population (default empty list; empty)
        """

        if group_id == "":
            vectors = self.populations[pop_id].vectors
        else:
            vectors = self.groups[group_id]

        self.populations[pop_id].addPathogensToVectors(genomes_numbers,vectors)

    def treatHosts(self, pop_id, frac_hosts, resistance_seqs, group_id=""):
        """Treat random fraction of infected hosts against some infection.

        Removes all infections with genotypes susceptible to given treatment.
        Pathogens are removed if they are missing at least one of the sequences
        in resistance_seqs from their genome. Removes this organism from
        population infected list and adds to healthy list if appropriate.

        Arguments:
        pop_id -- ID of population to be modified (String)
        frac_hosts -- fraction of hosts considered to be randomly selected
            (number between 0 and 1)
        resistance_seqs -- contains sequences required for treatment resistance
            (list of Strings)

        Keyword arguments:
        hosts -- list of specific hosts to sample from, if empty, samples from
            whole population (default empty list; empty)
        """

        if group_id == "":
            hosts = self.populations[pop_id].hosts
        else:
            hosts = self.groups[group_id]

        self.populations[pop_id].treatHosts(frac_hosts,resistance_seqs,hosts)

    def treatVectors(self, pop_id, frac_vectors, resistance_seqs, group_id=""):
        """Treat random fraction of infected vectors agains some infection.

        Removes all infections with genotypes susceptible to given treatment.
        Pathogens are removed if they are missing at least one of the sequences
        in resistance_seqs from their genome. Removes this organism from
        population infected list and adds to healthy list if appropriate.

        Arguments:
        pop_id -- ID of population to be modified (String)
        frac_vectors -- fraction of vectors considered to be randomly selected
            (number between 0 and 1)
        resistance_seqs -- contains sequences required for treatment resistance
            (list of Strings)

        Keyword arguments:
        vectors -- list of specific vectors to sample from, if empty, samples
            from whole population (default empty list; empty)
        """

        if group_id == "":
            vectors = self.populations[pop_id].vectors
        else:
            vectors = self.groups[group_id]

        self.populations[pop_id].treatVectors(
            frac_vectors,resistance_seqs,vectors
            )

    def protectHosts(
            self, pop_id, frac_hosts, protection_sequence, group_id=""):
        """Protect a random fraction of infected hosts against some infection.

        Adds protection sequence specified to a random fraction of the hosts
        specified. Does not cure them if they are already infected.

        Arguments:
        pop_id -- ID of population to be modified (String)
        frac_hosts -- fraction of hosts considered to be randomly selected
            (number between 0 and 1)
        protection_sequence -- sequence against which to protect (String)

        Keyword arguments:
        hosts -- list of specific hosts to sample from, if empty, samples from
            whole population (default empty list; empty)
        """

        if group_id == "":
            hosts = self.populations[pop_id].hosts
        else:
            hosts = self.groups[group_id]

        self.populations[pop_id].protectHosts(
            frac_hosts,protection_sequence,hosts
            )

    def protectVectors(
            self, pop_id, frac_vectors, protection_sequence, group_id=""):
        """Protect a random fraction of infected vectors against some infection.

        Adds protection sequence specified to a random fraction of the vectors
        specified. Does not cure them if they are already infected.

        Arguments:
        pop_id -- ID of population to be modified (String)
        frac_vectors -- fraction of vectors considered to be randomly selected
            (number between 0 and 1)
        protection_sequence -- sequence against which to protect (String)

        Keyword arguments:
        vectors -- list of specific vectors to sample from, if empty, samples
            from whole population (default empty list; empty)
        """

        if group_id == "":
            vectors = self.populations[pop_id].vectors
        else:
            vectors = self.groups[group_id]

        self.populations[pop_id].protectVectors(
            frac_vectors,protection_sequence,vectors
            )

    def setSetup(self, pop_id, setup_id):
        """Assign parameters stored in Setup object to this population.

        Arguments:
        pop_id -- ID of population to be modified (String)
        setup_id -- ID of setup to be assigned (String)
        """

        self.populations[pop_id].setSetup( self.setups[setup_id] )


    ### Preset fitness functions: ###

    @staticmethod
    def stabilizingSelection(genome, optimal_genome, min_fitness):
        """Evaluate genome fitness by decreasing with distance from optimal seq.

        A purifying selection fitness function based on exponential decay of
        fitness as genomes move away from the optimal sequence. Distance is
        measured as percent Hamming distance from an optimal genome sequence.

        Arguments:
        genome -- the genome to be evaluated (String)
        optimal_genome -- the genome sequence to measure distance against, has
            fitness of 1 (String)
        min_fitness -- minimum fitness value at maximum distance from optimal
            genome (number > 0)

        Return:
        fitness value of genome (number)
        """

        distance = td.hamming(genome, optimal_genome) / len(genome)
        fitness = np.exp( np.log( min_fitness ) * distance )

        return fitness

    @staticmethod
    def disruptiveSelection(genome, worst_genome, min_fitness):
        """Evaluate genome fitness by increasing with distance from worst seq.

        A purifying selection fitness function based on exponential decay of
        fitness as genomes move closer to the worst possible sequence. Distance
        is measured as percent Hamming distance from the worst possible genome
        sequence.

        Arguments:
        genome -- the genome to be evaluated (String)
        optimal_genome -- the genome sequence to measure distance against, has
            fitness of min_fitness (String)
        min_fitness -- fitness value of worst possible genome (number > 0)

        Return:
        fitness value of genome (number)
        """

        distance = td.hamming(genome, worst_genome) / len(genome)
        fitness = np.exp( np.log( min_fitness ) * ( 1 - distance ) )

        return fitness
