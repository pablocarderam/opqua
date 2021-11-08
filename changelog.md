
# Opqua Changelog

## v0.9.0
## 8 Nov 2021
All graphs in publication (title pending) generated with this stable version.

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
