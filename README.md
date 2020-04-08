Opqua
=====

**opqua** (opkua, upkua)
\[[Chibcha/muysccubun](https://en.wikipedia.org/wiki/Chibcha_language)\]

**I.** *noun*. ailment, disease, illness

**II.** *noun*. cause, reason \[*for which something occurs*\]

Taken from D. F. Gómez Aldana's
[muysca-spanish dictionary](http://muysca.cubun.org/opqua).

### Opqua is an epidemiological modeling framework for pathogen population genetics and evolution.

Opqua stochastically simulates pathogens with specific, evolving genotypes that
spread through populations of hosts that can have specific immune profiles.

Opqua is a useful tool to test out scenarios, explore hypotheses, and make
predictions about the relationship between pathogen evolution and epidemiology.

Among many things, Opqua can model
- host-host and vector-borne transmission,
- host recovery and death
- competition and evolution of pathogen strains across arbitrary adaptive
landscapes
- metapopulations with complex structure
- interventions altering demographic, ecological, or evolutionary parameters
- treatment and immunization of hosts or vectors

Opqua was developed by [Pablo Cárdenas](https://pablo-cardenas.com).
Follow my science antics at [@pcr_guy on Twitter](https://twitter.com/pcr_guy).

Opqua is available under an [MIT](https://choosealicense.com/licenses/mit/)
License.

## Installation


## Usage

All usage is handled through the Opqua ```Model``` class.

For example usage, have a look at the ```examples``` folder.

The Model class contains populations, setups, and interventions to be used
in simulation. Also contains groups of hosts/vectors for manipulations and
stores model history as snapshots for each time point.

## Attributes
populations -- dictionary with keys=population IDs, values=Population
    objects
setups -- dictionary with keys=setup IDs, values=Setup objects
interventions -- contains model interventions in the order they will occur
groups -- dictionary with keys=group IDs, values=lists of hosts/vectors
self.history -- dictionary with keys=time values, values=Model objects that
    are snapshots of Model at that timepoint

## Methods list

### Model Initialization and Simulation: ###

- newSetup -- creates a new Setup, save it in setups dict under given name
- newIntervention -- creates a new intervention executed during simulation
- run -- simulates model for a specified length of time

### Data Output and Plotting: ###

- saveToDataFrame -- saves status of model to dataframe, writes to file
- getPathogens -- creates Dataframe with counts for all pathogen genomes
- getProtections -- creates Dataframe with counts for all protection sequences
- populationsPlot -- plots aggregated totals per population across time
- compartmentPlot -- plots number of naive, infected, recovered, dead
hosts/vectors vs time
- compositionPlot -- plots counts for pathogen genomes or resistance vs. time

### Model interventions: ###

#### Make and connect populations ####
- newPopulation -- create a new Population object with setup parameters
- linkPopulations -- set migration rate from one population towards another
- createInterconnectedPopulations -- create new populations, link all of them
    to each other

#### Modify population parameters ####
- setSetup -- assigns a given set of parameters to this population

#### Manipulate hosts and vectors in population ####
- newHostGroup -- returns a list of random (healthy or any) hosts
- newVectorGroup -- returns a list of random (healthy or any) vectors
- addHosts -- adds hosts to the population
- addVectors -- adds vectors to the population
- removeHosts -- removes hosts from the population
- removeVectors -- removes vectors from the population
- addPathogensToHosts -- adds pathogens with specified genomes to hosts
- addPathogensToVectors -- adds pathogens with specified genomes to vectors
- treatHosts -- removes infections susceptible to given treatment from hosts
- treatVectors -- removes infections susceptible to treatment from vectors
- protectHosts -- adds protection sequence to hosts
- protectVectors -- adds protection sequence to vectors


### Preset fitness functions: ###

- stabilizingSelection -- evaluates genome fitness by decreasing with distance
    from optimal sequence
- disruptiveSelection -- evaluates genome fitness by increasing with distance
    from worst sequence

## Detailed list:

```Model()```

Create a new Model object.

```newSetup()```

Create a new Setup, save it in setups dict under given name.

Two preset setups exist: "vector-borne" and "host-host". You may select
one of the preset setups with the preset keyword argument and then
modify individual parameters with additional keyword arguments, without
having to specify all of them.

Preset parameter setups:

"vector-borne":
  num_loci = 10
  possible_alleles = 'ATCG'
  fitnessHost = (lambda g: 1)
  fitnessVector = (lambda g: 1)
  contact_rate_host_vector = 1e1
  contact_rate_host_host = 0
  mean_inoculum_host = 1e2
  mean_inoculum_vector = 1e2
  recovery_rate_host = 1e-1
  recovery_rate_vector = 1e-2
  recombine_in_host = 0
  recombine_in_vector = 1e-2
  mutate_in_host = 1e-6
  mutate_in_vector = 0
  death_rate_host = 0
  death_rate_vector = 0
  protection_upon_recovery_host = None
  protection_upon_recovery_vector = None

"host-host":
  num_loci = 10
  possible_alleles = 'ATCG'
  fitnessHost = (lambda g: 1)
  fitnessVector = (lambda g: 1)
  contact_rate_host_vector = 0
  contact_rate_host_host = 2e1
  mean_inoculum_host = 1e1
  mean_inoculum_vector = 0
  recovery_rate_host = 1e-1
  recovery_rate_vector = 1e1
  recombine_in_host = 1e-3
  recombine_in_vector = 0
  mutate_in_host = 1e-6
  mutate_in_vector = 0
  death_rate_host = 0
  death_rate_vector = 0
  protection_upon_recovery_host = None
  protection_upon_recovery_vector = None

Arguments:
- name -- name of setup to be used as a key in model setups dictionary

Keyword arguments:
- preset -- preset setup to be used: "vector-borne" or "host-host", if
None, must define all other keyword arguments (default None; None or
String)
- num_loci -- length of each pathogen genome string (int > 0)
possible_alleles -- set of possible characters in all genome string, or
at each position in genome string (String or list of Strings with
- num_loci elements)
- fitnessHost -- relative fitness in head-to-head competition within host
(number >= 0)
- fitnessVector -- relative fitness in head-to-head competition within
vector (number >= 0)
- contact_rate_host_vector -- rate of host-vector contact events, not
necessarily transmission, assumes constant population density;
evts/time (number >= 0)
- contact_rate_host_host -- rate of host-host contact events, not
necessarily transmission, assumes constant population density;
evts/time (number >= 0)
- mean_inoculum_host -- mean number of pathogens that are transmitted from
a vector or host into a new host during a contact event (int >= 0)
- mean_inoculum_vector -- mean number of pathogens that are transmitted
from a host to a vector during a contact event (int >= 0)
- recovery_rate_host -- rate at which hosts clear all pathogens;
1/time (number >= 0)
- recovery_rate_vector -- rate at which vectors clear all pathogens
1/time (number >= 0)
- recombine_in_host -- rate at which recombination occurs in host;
evts/time (number >= 0)
- recombine_in_vector -- rate at which recombination occurs in vector;
evts/time (number >= 0)
- mutate_in_host -- rate at which mutation occurs in host; evts/time
(number >= 0)
- mutate_in_vector -- rate at which mutation occurs in vector; evts/time
(number >= 0)
- death_rate_host -- infected host death rate; 1/time (number >= 0)
- death_rate_vector -- infected vector death rate; 1/time (number >= 0)
- protection_upon_recovery_host -- defines indexes in genome string that
define substring to be added to host protection sequences after
recovery (None or array-like of length 2 with int 0-num_loci)
- protection_upon_recovery_vector -- defines indexes in genome string that
define substring to be added to vector protection sequences after
recovery (None or array-like of length 2 with int 0-num_loci)
