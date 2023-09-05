# `Model` Documentation

All usage is handled through the Opqua `Model` class.
The `Model` class contains populations, setups, and interventions to be used
in simulation. It also contains groups of hosts/vectors for manipulations and
stores model history as snapshots for specific time points.

To use it, import the class as

```python
from opqua.model import Model
```

You can find a detailed account of everything `Model` does in the
[Model attributes](#model-class-attributes) and
[Model class methods list](#model-class-methods-list) sections.

## `Model` class attributes

- `populations` -- dictionary with keys=population IDs, values=Population
    objects
- `setups` -- dictionary with keys=setup IDs, values=Setup objects
- `interventions` -- contains model interventions in the order they will occur
- `groups` -- dictionary with keys=group IDs, values=lists of hosts/vectors
- `history` -- dictionary with keys=time values, values=Model objects that
    are snapshots of Model at that timepoint
- `global_trackers` -- dictionary keeping track of some global indicators over all
    the course of the simulation
- `custom_condition_trackers` -- dictionary with keys=ID of custom condition,
    values=functions that take a Model object as argument and return True or
    False; every time True is returned by a function in
    custom_condition_trackers, the simulation time will be stored under the
    corresponding ID inside global_trackers['custom_condition']
- `t_var` -- variable that tracks time in simulations

The dictionary global_trackers contains the following keys:
- `num_events`: dictionary with the number of each kind of event in the simulation
- `last_event_time`: time point at which the last event in the simulation happened
- `genomes_seen**: list of all unique genomes that have appeared in the
    simulation
- `custom_conditions`: dictionary with keys=ID of custom condition, values=lists
    of times; every time True is returned by a function in
    custom_condition_trackers, the simulation time will be stored under the
    corresponding ID inside global_trackers['custom_condition']

The dictionary `num_events` inside of global_trackers contains the following keys:
- `MIGRATE_HOST`
- `MIGRATE_VECTOR`
- `POPULATION_CONTACT_HOST_HOST`
- `POPULATION_CONTACT_HOST_VECTOR`
- `POPULATION_CONTACT_VECTOR_HOST`
- `CONTACT_HOST_HOST`
- `CONTACT_HOST_VECTOR`
- `CONTACT_VECTOR_HOST`
- `RECOVER_HOST`
- `RECOVER_VECTOR`
- `MUTATE_HOST`
- `MUTATE_VECTOR`
- `RECOMBINE_HOST`
- `RECOMBINE_VECTOR`
- `KILL_HOST`
- `KILL_VECTOR`
- `DIE_HOST`
- `DIE_VECTOR`
- `BIRTH_HOST`
- `BIRTH_VECTOR`

KILL_HOST and KILL_VECTOR denote death due to infection, whereas DIE_HOST and
DIE_VECTOR denote death by natural means.

## `Model` class methods list

### Model initialization and simulation

- `setRandomSeed()` -- set random seed for numpy random number
generator
- `newSetup()` -- creates a new Setup, save it in setups dict under
given name
- `newIntervention()` -- creates a new intervention executed
during simulation
- `run()` -- simulates model for a specified length of time
- `runReplicates]()` -- simulate replicates of a model, save only
end results
- `runParamSweep()` -- simulate parameter sweep with a model, save
only end results
- `copyState()` -- copies a slimmed-down representation of model state
- `deepCopy()` -- copies current model with inner references

### Data Output and Plotting

- `saveToDataFrame()` -- saves status of model to data frame,
writes to file
- `getPathogens()` -- creates data frame with counts for all
pathogen genomes
- `getProtections()` -- creates data frame with counts for all
protection sequences
- `populationsPlot()` -- plots aggregated totals per
population across time
- `compartmentPlot()` -- plots number of naive, infected,
recovered, dead hosts/vectors vs time
- `compositionPlot()` -- plots counts for pathogen genomes or
resistance vs. time
- `clustermap()` -- plots heatmap and dendrogram of all pathogens in
given data
- `pathogenDistanceHistory()` -- calculates pairwise
distances for pathogen genomes at different times
- `getGenomeTimes()` -- create DataFrame with times genomes first
appeared during simulation
- `getCompositionData()` -- create dataframe with counts for
      pathogen genomes or resistance

### Model interventions

#### Make and connect populations:
- `newPopulation()` -- create a new Population object with
setup parameters
- `linkPopulationsHostMigration()` -- set host
migration rate from one population towards another
- `linkPopulationsVectorMigration()` -- set
vector migration rate from one population towards another
- `linkPopulationsHostHostContact()` -- set
host-host inter-population contact rate from one population towards another
- `linkPopulationsHostVectorContact()` -- set
host-vector inter-population contact rate from one population towards another
- `linkPopulationsVectorHostContact()` -- set
vector-host inter-population contact rate from one population towards another
- `createInterconnectedPopulations()` -- create new populations, link all of them to 
each other by migration and/or inter-population contact

#### Manipulate hosts and vectors in population:
- `newHostGroup()` -- returns a list of random (healthy or any)
hosts
- `newVectorGroup()` -- returns a list of random (healthy or
  any) vectors
- `addHosts()` -- adds hosts to the population
- `addVectors()` -- adds vectors to the population
- `removeHosts](#removehosts)` -- removes hosts from the population
- `removeVectors()` -- removes vectors from the population
- `addPathogensToHosts()` -- adds pathogens with
specified genomes to hosts
- `addPathogensToVectors()` -- adds pathogens with
specified genomes to vectors
- `treatHosts()` -- removes infections susceptible to given
treatment from hosts
- `treatVectors()` -- removes infections susceptible to
treatment from vectors
- `protectHosts()` -- adds protection sequence to hosts
- `protectVectors()` -- adds protection sequence to vectors
- `wipeProtectionHosts()` -- removes all protection
sequences from hosts
- `wipeProtectionVectors()` -- removes all protection
sequences from vectors

#### Modify population parameters:
- `setSetup()` -- assigns a given set of parameters to this population

#### Utility:
- `customModelFunction()` -- returns output of given function run on model

### Preset fitness functions

- `peakLandscape()` -- evaluates genome numeric phenotype by
decreasing with distance from optimal sequence
- `valleyLandscape()` -- evaluates genome numeric phenotype by
increasing with distance from worst sequence


## Detailed `Model` documentation

```{eval-rst}
.. autoclass:: opqua.model.Model
   :members:
```