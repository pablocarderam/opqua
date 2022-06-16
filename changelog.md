
# Opqua Changelog

## v0.9.8-immunity-0.6
## 16 Jun 2022
Fixed two bugs related to excessive immunity loss and added print_every_n_events
parameter to run() function.

- In the first one, host/vector IMMUNE coefficients were set to zero upon recovery
regardless of actual immune state. Fixed by moving healthyCoefficientRow() to
Host and Vector classes, then modifying the function to check that host/vector's
immune state and set the IMMUNE coefficient accordingly.

- In the second bug, loss of immunity happened upon migration because immune
sequences were not being transferred upon migration. Fixed that in migrate().

- Also added print_every_n_events parameter to run() function in Model class.

## v0.9.8-immunity-0.5
## 17 May 2022
Change immunization architecture slightly. Now, coefficient matrix has an
additional "IMMUNIZED" column tracking whether hosts/vectors have any immunity
at all. This is used to identify hosts/vectors capable of being deimmunized.
Additionally, individual coefficients are now recalculated in all instances
where immunity changes.

## v0.9.8-immunity-0.4
## 20 Apr 2022
Change all mentions of "lethality" to "mortality".

## v0.9.8-immunity-0.3
## 13 Mar 2022
Update Pillow to 9.0.1

## v0.9.8-immunity-0.2
## 9 Mar 2022
Major overhaul of immunity mechanics in simulation engine. Still testing some
behaviors, but going well so far.

- Added new events for acquisition and loss of immunity; immunity is no longer
  acquired only upon recovery.
- Added new parameters to control acquisition and loss of immunity.
- Added new parameter to specify custom user functions to evaluate effect of
  immunity by comparing pathogen genomes with genomes in immune memory.
- Immune memory now samples entire genome (so "partial matches" can be specified
  through custom functions as explained in the last point).
- Immunity now consists of modifying intra-host fitness through an extra
  coefficient; if it is "sterilizing" (i.e. coefficient=0) then that pathogen is
  removed from host/vector.
- Immunity therefore affects rates of all event types through a pathogen's
  intra-host fitness, but it can additionally decrease rates of events where
  pathogens are doing the action (transmission, death, mutation, recombination,
  immunization)
- Immunity can be inherited vertically on a sequence-by-sequence basis, not
  all-or-nothing
- Added an "immunity" example to the examples/tutorials folder.

Also added other modifications done in parallel on the main branch:
- Modify behavior of compositionPlot to take all hosts/vectors into account when
  plotting immunity/protection.

- Change units of lethality_rate_host and lethality_rate_vector to be a rate
  like other parameters for better internal consistency, instead of a fraction
  of recovered cases. Affects only behavior of simulations with disease-caused
  mortality. Updated example accordingly.

- Fix bug in updateVectorCoefficients() specific to natality and migration.
  Does not affect any simulation results as long as natality and migration are
  not functions of pathogen genome sequence.

## v0.9.8-immunity-0.1
## 2 Mar 2022
Renamed "Protection" with "Immunity".

## v0.9.8
## 27 Feb 2022
Update compositionPlot additional arguments.

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
infinite loops when tampering with t_var.
- changed order in which time delta was added to t_var to be after interventions
  occur
- do not carry out interventions if t_var is past simulation end time

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
- made time variable t_var a property of Model objects instead of internal to
  Gillespie object, allowing users to modify simulation time (e.g. for
  killswitches); does not affect simulation results

Opqua structure changes:
- changed way interventions are executed to guarantee that the Model object
  being simulated carries the intervention upon itself
- changed runReplicates() to create independent copies of model object and run
  simulations on each (should not affect simulation results due to
  parallelization through joblib)
- added deepCopy() to reassign all internal model and population references in
  copied model objects
- changed runReplicates() and runParamSweep() to use deepCopy() (should not
  affect simulation results due to parallelization through joblib)
- added customModelFunction() function to allow users to add custom methods to
  specific Model instances (e.g. for killswitches and conditional interventions)

Opqua syntax changes:
- changed syntax of interventions to force all functions used to be methods of
  the Model object being intervened

## v0.9.0
## 8 Nov 2021
All graphs in publication (title pending) generated with this stable version.
[RETROACTIVE EDIT: first draft only]

General model structure changes:
- added transmission_efficiency_host_host, transmission_efficiency_host_vector,
  transmission_efficiency_vector_host as additional parameters
- made global_trackers copy into history object of a model
- adjusted computation of recombination probabilities for hosts and vectors
- skip recombining when parental genomes are the same
- changed genome sampling during inoculation of hosts and vectors
- added more flexibility, changed structure, and debugged migration/population
  contact options in runParamSweep()

In compositionDf():
- remove missing data
- change algorithm to replace combinations of genomes (more efficient when
  combinations are limiting factor)
- allow user to specify genome plotting order

Miscellaneous:
- add option to plot population fractions instead of absolute counts in
  compositionPlot()
- Gillespie algorithm now prints out event name rather than ID number
- changed error handling when adding pathogens to hosts and vectors
- removed a duplicate definition in Model newSetup()
- changed a group name in intervention_example.py

## v0.2.6
## 8 Oct 2021
v0.2.5 created a major bug that escaped my attention with the division by zero
error fix.

- corrected Host and Vector acquirePathogen() functions to restore correct
      behavior
- added a requirements.txt file purely for reference purposes in case a future
      dependency update breaks opqua

## v0.2.5
## 8 Oct 2021
- added global_trackers dictionary to Model in order to track and return some
      global indicators when running replicates or parameter sweeps
- added addCustomConditionTracker() function to Model class in order to allow
      users to track custom events in model
- modified mutation() and recombination() in Vector and Host classes as well as
      addPathogensToHosts() and addPathogensToVectors() in Population class in
      order to track genomes seed for global_trackers
- added model attribute to Population class for the above reason as well
- fix way host and vector sampling is handled (apparently I forgot to actually
      implement it in Gillespie and Population classes after I added the
      arguments, haha)
- minor bug fix: modify acquirePathogen in Vector and Host classes to avoid
      division by zero errors when recalculating sum_fitness after it was zero
- corrected name of compositionDataframe argument in compositionPlot() to
      composition_dataframe

## v0.2.4
## 5 Oct 2021
- fixed regex processing bug in compositionDf()
- added **kwargs argument passing to joblib functions to allow user to change
      backend and stuff

Trying to deploy on cluster so bear with me on the updates here

## v0.2.3
## 5 Oct 2021
- added parameter sweep function runParamSweep()
- added id property and argument to Setup() in order to associate a Setup to its
      ID in a Model, so that runParamSweep() can edit the setups
- added getCompositionData() function to Model class to allow user output
      composition data without plotting compartments
- fixed bug in how runReplicates() computed and return output
- added verbose optional argument to saveToDf() to reduce console output
- added composition_dataframe optional argument to allow for pre-computed data
- added setup.cfg as per [Joel Barmettler](https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56)
      I guess?

## v0.2.2
## 21 July 2021
- change compositionPLot remove_legend behavior to fix bug
- change pathogenDistanceDf seq_names behavior to fix bug
- reduce mean inoculum from hosts into vectors to reflect malaria cycle
- modify infectHost and infectVector inoculation behavior so that mean_inoculum
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
- make RECEIVE_CONTACT and RECEIVE_POPULATION_CONTACT coefficient columns
      in arrays, modify how these events are handled. Done, tested
- migrate vectors. Done, test
- contact between populations (without migration). Done, tested

### Output
- make genome labels optional, if no labels write a file with the genomes
      in the same order for compositionPlot. Done, test
- compositionPlot – make custom groupings, eg. "genomes containing
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
- correctly update arguments in function documentation and README for
compositionPlot family functions
