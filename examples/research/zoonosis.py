
"""
Models SARS-CoV-2 emergence.
"""

from opqua.model import Model

bat_genome = 'GENOMEBAT_BATSPIKE' # optimal bat genome
pgn_genome = 'GENOMEPGN_PGNSPIKE' # optimal pangolin genome
hmn_genome = 'GENOMEBAT_PGNSPIKE' # optimal human genome

def batFitness(genome):
    return Model.stabilizingSelection(
        genome, optimal_genome=bat_genome, min_fitness=1e-10
        )

def pgnFitness(genome):
    return Model.stabilizingSelection(
        genome, optimal_genome=pgn_genome, min_fitness=1e-10
        )

def batFitness(genome):
    return Model.stabilizingSelection(
        genome, optimal_genome=hmn_genome, min_fitness=1e-10
        )

m = Model()
m.newSetup( # Now, we'll define our new setup:
    'setup', preset='host-host', # Use default host-host parameters.
    possible_alleles=set(bat_genome+pgn_genome+hmn_genome),
        # Define "letters" in the "genome", or possible alleles for each locus.
    num_loci=len(bat_genome),
        # Define length of "genome", or total number of alleles.
    mutate_in_host=0, # Modify mutation rate
    recombine_in_host=2e-2 # Modify mutation rate
    )

bats_pop = 1e2
pgns_pop = 1e2
hmns_pop = 1e2
ref_pop = 1e2
std_mig_rate = 1e-3

m.newPopulation('Bats','setup',num_hosts=bats_pop)
m.newPopulation('Pangolins','setup',num_hosts=pgns_pop)
m.newPopulation('Humans','setup',num_hosts=hmns_pop)

m.linkPopulations('Bats','Pangolins',std_mig_rate*ref_pop/bats_pop)
m.linkPopulations('Pangolins','Bats',std_mig_rate*ref_pop/pgns_pop)
m.linkPopulations('Bats','Humans',std_mig_rate*ref_pop/bats_pop)
m.linkPopulations('Humans','Bats',std_mig_rate*ref_pop/hmns_pop)
m.linkPopulations('Pangolins','Humans',std_mig_rate*ref_pop/pgns_pop)
m.linkPopulations('Humans','Pangolins',std_mig_rate*ref_pop/hmns_pop)

m.addPathogensToHosts( 'Bats',{bat_genome:bats_pop/10} )
m.addPathogensToHosts( 'Pangolins',{pgn_genome:pgns_pop/10} )

m.protectHosts('Humans',1,'T_B')
m.protectHosts('Humans',1,'N_P')

m.run(0,500)
data = m.saveToDataFrame('zoonosis.csv')

zoonosis_composition_human = m.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'zoonosis_composition_human.png', data,populations=['Humans'],
    num_top_sequences=7
    )
zoonosis_composition_bat = m.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'zoonosis_composition_bat.png', data,populations=['Bats'],
    num_top_sequences=7
    )
zoonosis_composition_pgn = m.compositionPlot(
        # Create a plot to track pathogen genotypes across time.
    'zoonosis_composition_pangolin.png', data,populations=['Pangolins'],
    num_top_sequences=7
    )

graph = m.populationsPlot( # Plot infected hosts per population over time.
    'zoonosis_population.png', data,
    y_label='Infected hosts' # change y label
    )
