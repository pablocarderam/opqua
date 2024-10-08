
# Opqua Changelog

## 23 Mar 2024
Major overhaul of computation of establishment and emergence probabilities in
`Landscape` class. This new approach is based on discretizing mutant cloud
diffusion to compute emergence probability at different pseudotimes (mutant
cloud sizes) after bottleneck, then computing probability of fixation within
pseudotime (mutant cloud size) frame, then normalizing by time to obtain true
rate (slideshow explains, preprint will explain in detail). Now each mutant has
a series of emergence-establishment rates corresponding to different mutant
cloud sizes/pseudotimes after bottleneck. Plot class changed accordingly to
visualize sum of all rates across time.

## 11 Feb 2024
Changed computation of emergence probabilities in `Landscape` class to correctly
account for pathway dependence in instances in which intermediates with lower
fitness are present. Also corrected wait time before mutations are accumulated.
Updated parameters, documentation, tutorial, and output in corresponding manner.

## 5 Feb 2024
Bug fixes and improvements for landscape mapping:
- incomplete mapping due to incorrect depth handling in `evaluateNeighbors()`
- incorrect calculation of probabilities of acquiring specific mutations: first,
  the total number of alleles used in normalization was wrong; second, the total
  number of available paths to acquiring a group of mutations (calculated by
  permutation without replacement and *assuming equal path likelihood*) was
  missing
- added fitness as a recorded variable in mutation networks and visualization
- improved visualization controls
- cleaned up documentation of `Landscape` class and `visualizeMutationNetwork`
- added landscape mapping example to `examples` folder and output graph to `img`
- split some of the landscaping parameters into host and vector versions

Also, added all necessary parameters for new intrahost evolution algorithm.

## 3 Feb 2024
Successfully implemented Landscape class to traverse and compute fitness
landscapes, the parameters needed for this in the `Setup` class, as well as a
plotting function to visualize the resulting mutation networks using the PyVis
package.

## v1.2.1
## 1 Feb 2024
Fixed some bugs in function importing (or lack thereof; recommended practice is 
to declare function parameters as Python code in same scope as simulation and
pass function names as arguments to `newSetup()` or `loadSetup()`).

Also corrected small bug in which `importlib_resources` wasn't being imported in
`Setup` class.

## v1.2.0
## 1 Feb 2024
Woke up and decided yesterday's changes were significant enough to warrant a
version change to 1.2.0 :)

## v1.1.2
## 1 Feb 2024
Forgot the `MANIFEST.in` file to get the csv files.

## v1.1.1
## 31 Jan 2024
- Added methods to save and load parameter settings into `Setup` objects using
  CSV files
- Changed `Setup` object internal mechanics
- Added default parameter values as external csv files, removed them from `Model`
  class

## 18 Jan 2024
Moved plot and data files to a separate directory for organizational purposes.

## v1.1.0
## 17 Jan 2024
Update to pandas dataframe handling (replace append with concat) to comply with
new version.

Cleaned up everything for a release. Looking ahead:
1. A big overhaul of mutation mechanisms and intra-host dynamics
2. A big overhaul of host acquired immunity
3. Performance improvements by refactoring simulation engine in C

## May 2023
Changes:
- Uses new Numpy Generator class for random numbers, some performance
  improvement expected.
- Added function ``zeroTruncatedPoisson`` to return random numbers from a
  zero-truncated Poisson distribution, to be used to define the number of
  cross-over events
- Added function ``zeroTruncatedNegBinomial`` to return random numbers from a
  zero-truncated negative binomial distribution (or at least an approximation),
  to be used to define inoculum size during transmission (Sobel et al. 2017 find
  negative binomial is superior to Poisson when estimating bottleneck size);
  also add model parameters ``variance_inoculum_host`` and
  ``variance_inoculum_vector`` to be used by this function
- Renamed `Gillespie` class as `Simulation`; added new simulation algorithm that
  implements a variation of tau leaping (but left exact Gillespie as default
  since estimating adequate tau parameter seems tricky)
- Fixed bug in `compartmentDf` which created duplicated rows in dataframe
  when hosts or vectors died infected and/or protected

## 3 Jan 2023
- Bumped Pillow to satisfy the Github bot
- Fixed a bug that made recombination events depend on mutation coefficients
  instead of recombination coefficients (does not affect published results)
- Slightly alter the way Poisson distributions are used to define the number of
  cross-over events and pathogens inoculated in transmission, for consistency
  (add 1 to the mean of events once they are guaranteed to happen; impact on
  existing simulations is negligible)
- Modify ``getWeightedRandom()`` to use numpy arrays (profiling shows it does
  increase efficiency)
- Added a new parameter to the ``run()`` function: ``skip_uninfected``. When
  ``True``, allows Opqua to store copies of only infected hosts/vectors as
  simulation progresses, and stores total number of healthy hosts to then
  reconstitute those as generic rows on dataframe after simulation is over. For
  simulations with a large number of hosts/vectors and relatively low infection
  prevalence, this can greatly increase simulation speed.
- Changed Gillespie algorithm to only recalculate probabilities for events that
  may happen in simulation (since many simulations omit certain types of events)
- Added parameters to `setSetup` that optionally allow recalculation of all
  host and/or vector coefficients, thus overwriting all establishment frequency
  effects; the new default is to not recalculate

## v1.0.2
## 17 Nov 2022
Previous fix didn't cut it and I jumped the gun on the release. This works.

## v1.0.1
## 17 Nov 2022
Fixed recombination bug where one of the two progeny genomes was being lost
(thanks David Suárez!).
Updated joblib version.

## v1.0.0
## 10 May 2022
Revisions incorporated!

## v0.9.9
## 20 Apr 2022
Change all mentions of "lethality" to "mortality".

## v0.9.8.4
## 13 Mar 2022
Update Pillow to 9.0.1

## v0.9.8.3
## 9 Mar 2022
Modify behavior of compositionPlot to take all hosts/vectors into account when
plotting immunity/protection.

## v0.9.8.2
## 8 Mar 2022
Change units of `lethality_rate_host` and `lethality_rate_vector` to be a rate like
other parameters for better internal consistency, instead of a fraction of
recovered cases. Affects only behavior of simulations with disease-caused
mortality. Updated example accordingly.

## v0.9.8.1
## 7 Mar 2022
Fix bug in `updateVectorCoefficients()` specific to natality and migration.
Does not affect any simulation results as long as natality and migration are not
functions of pathogen genome sequence.

## v0.9.8
## 27 Feb 2022
Update `compositionPlot` additional arguments.

## v0.9.7
## 18 Jan 2022
Update Pillow to 9.0.0.

## 21 Dec 2021
Added biorXiv links.

## v0.9.6
## 14 Dec 2021
Changed computation of inter-population contact rates to match logic of
intra-population contact rates. Does not affect outcome of simulations not using
inter-population contact.

## v0.9.5
## 9 Dec 2021
Parenthesis error.

## v0.9.4
## 9 Dec 2021
Another small change to the Gillespie algorithm, this time to avoid rare
infinite loops when tampering with `t_var`.
- changed order in which time delta was added to `t_var` to be after interventions
  occur
- do not carry out interventions if `t_var` is past simulation end time

## v0.9.3
## 8 Dec 2021
Same bug as below, the fix was incomplete.

## v0.9.2
## 8 Dec 2021
There was a bug in the last release regarding interventions! in order to
correctly impplement custom user killswitches, it is important to update the
time variable t_var immediately *before* an intervention takes place, not after
it. Fixed in Gillespie method now, does not change behavior of any previous
simulations.

## v0.9.1
## 8 Dec 2021

General simulation structure changes:
- changed handling of intra-population contact rate between hosts and vectors
  within the Gillespie function to fit a biting rate definition for the contact
  rate (i.e. constant per vector contact rate). Does not affect behavior of
  simulations if total number of hosts and vectors is equal (or if models have
  no vector-borne pathogens)
- made time variable t_var a property of `Model` objects instead of internal to
  `Gillespie` object, allowing users to modify simulation time (e.g. for
  killswitches); does not affect simulation results

Opqua structure changes:
- changed way interventions are executed to guarantee that the Model object
  being simulated carries the intervention upon itself
- changed `runReplicates()` to create independent copies of model object and run
  simulations on each (should not affect simulation results due to
  parallelization through joblib)
- added `deepCopy()` to reassign all internal model and population references in
  copied model objects
- changed `runReplicates()` and `runParamSweep()` to use `deepCopy()` (should not
  affect simulation results due to parallelization through joblib)
- added `customModelFunction()` function to allow users to add custom methods to
  specific Model instances (e.g. for killswitches and conditional interventions)

Opqua syntax changes:
- changed syntax of interventions to force all functions used to be methods of
  the Model object being intervened

## v0.9.0
## 8 Nov 2021
All graphs in publication (title pending) generated with this stable version.
[RETROACTIVE EDIT: first draft only]

General model structure changes:
- added transmission_efficiency_host_host, `transmission_efficiency_host_vector`,
  `transmission_efficiency_vector_host` as additional parameters
- made global_trackers copy into history object of a model
- adjusted computation of recombination probabilities for hosts and vectors
- skip recombining when parental genomes are the same
- changed genome sampling during inoculation of hosts and vectors
- added more flexibility, changed structure, and debugged migration/population
  contact options in `runParamSweep()`

In compositionDf():
- remove missing data
- change algorithm to replace combinations of genomes (more efficient when
  combinations are limiting factor)
- allow user to specify genome plotting order

Miscellaneous:
- add option to plot population fractions instead of absolute counts in
  `compositionPlot()`
- Gillespie algorithm now prints out event name rather than ID number
- changed error handling when adding pathogens to hosts and vectors
- removed a duplicate definition in `Model newSetup()`
- changed a group name in `intervention_example.py`

## v0.2.6
## 8 Oct 2021
v0.2.5 created a major bug that escaped my attention with the division by zero
error fix.

- corrected `Host` and `Vector` `acquirePathogen()` functions to restore correct
      behavior
- added a requirements.txt file purely for reference purposes in case a future
      dependency update breaks opqua

## v0.2.5
## 8 Oct 2021
- added global_trackers dictionary to Model in order to track and return some
      global indicators when running replicates or parameter sweeps
- added addCustomConditionTracker() function to `Model` class in order to allow
      users to track custom events in model
- modified `mutation()` and `recombination()` in `Vector` and `Host` classes as well as
      `addPathogensToHosts()` and `addPathogensToVectors()` in `Population` class in
      order to track genomes seed for `global_trackers`
- added model attribute to `Population` class for the above reason as well
- fix way host and vector sampling is handled (apparently I forgot to actually
      implement it in `Gillespie` and `Population` classes after I added the
      arguments, haha)
- minor bug fix: modify `acquirePathogen` in `Vector` and `Host` classes to avoid
      division by zero errors when recalculating `sum_fitness` after it was zero
- corrected name of `compositionDataframe` argument in `compositionPlot()` to
      composition_dataframe

## v0.2.4
## 5 Oct 2021
- fixed regex processing bug in `compositionDf()`
- added `**kwargs` argument passing to joblib functions to allow user to change
      backend and stuff

Trying to deploy on cluster so bear with me on the updates here

## v0.2.3
## 5 Oct 2021
- added parameter sweep function `runParamSweep()`
- added id property and argument to `Setup()` in order to associate a Setup to its
      ID in a Model, so that `runParamSweep()` can edit the setups
- added `getCompositionData()` function to `Model` class to allow user output
      composition data without plotting compartments
- fixed bug in how `runReplicates()` computed and return output
- added verbose optional argument to `saveToDf()` to reduce console output
- added composition_dataframe optional argument to allow for pre-computed data
- added setup.cfg as per [Joel Barmettler](https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56)
      I guess?

## v0.2.2
## 21 July 2021
- change `compositionPLot` `remove_legend` behavior to fix bug
- change `pathogenDistanceDf` `seq_names` behavior to fix bug
- reduce mean inoculum from hosts into vectors to reflect malaria cycle
- modify `infectHost` and `infectVector` inoculation behavior so that mean_inoculum
      does not affect overall transmission rate; each infection now results in
      at least 1 pathogen transfer (if not containing and not immune to the
      pathogen genome sampled)

## v0.2.1
## 1 June 2021
- update version tags

## v0.2.0
## 1 June 2021
### Simulation mechanics
- pathogen genome influences transmission, death, recovery, migration,
      mutation probabilities. Done, tested
- independent recombination of alleles -> make chromosome separators!!
      make reassortment parameter option as related but separate to
      recombination. Done, tested
- host/vector birth/death rates in populations -> make birth event.
      Done, tested
- separate natural death into an event that doesn't log deaths. Done, tested
- make `RECEIVE_CONTACT` and `RECEIVE_POPULATION_CONTACT` coefficient columns
      in arrays, modify how these events are handled. Done, tested
- migrate vectors. Done, test
- contact between populations (without migration). Done, tested

### Output
- make genome labels optional, if no labels write a file with the genomes
      in the same order for `compositionPlot`. Done, test
- `compositionPlot` – make custom groupings, eg. "genomes containing
      sequence AAA". Done, test
- make option to count only 1 fitness-dominant strain/host in
      compositionPlot. Done, test
- genomes and dates output for TDA. Done, test

### Technical details
- set seed. Done, test
- make all events except mutation, recombination and recovery
      1-coefficient like contact and pop contact are. this way, coefficients are
      fraction of healthy rate. Modify dummy rows in coeff arrays to match
      this. Done, tested
- parallelizeable simulations. Done, test
- try Numba, Cython, JAX optimization

### Miscellaneous
- update docs. Done
- update tutorials. Done, tested

### TODO:
- correctly update arguments in function documentation and `README` for
`compositionPlot` family functions
