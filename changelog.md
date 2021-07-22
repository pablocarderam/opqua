
# Opqua Changelog

## v0.2.2
## 21 July 2021
- change compositionPLot remove_legend behavior to fix bug
- change pathogenDistanceDf seq_names behavior to fix bug
- reduce mean inoculum from hosts into vectors to reflect malaria cycle
- modify infectHost and infectVector inoculation behavior so that mean_inoculum
does not affect overall transmission rate; each infection now results in at
least 1 pathogen transfer (if not containing and not immune to the pathogen
genome sampled)

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
