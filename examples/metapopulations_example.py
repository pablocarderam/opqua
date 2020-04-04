
'''
Simple model
'''

from opqua.model import Model

model = Model()
model.newSetup('setup_normal', default='vector-borne') # uses default parameters
model.newSetup('setup_cluster', contact_rate_host_vector=2e1 ,default='vector-borne') # uses default parameters

model.newPopulation('population_A','setup_normal', num_hosts=20, num_vectors=20)
model.newPopulation('population_B','setup_normal', num_hosts=20, num_vectors=20)
model.newPopulation('isolated_population','setup_normal', num_hosts=20, num_vectors=20)

model.createInterconnectedPopulations(5,1,'clustered_population_','setup_cluster', num_hosts=20, num_vectors=20)
model.linkPopulations('population_A','clustered_population_4',1)
model.linkPopulations('population_A','population_B',1)

model.addPathogensToHosts( 'population_A',{'AAAAAAAAAA':5} )
model.addPathogensToHosts( 'population_B',{'GGGGGGGGGG':5} )

output = model.run(0,100)
data = model.saveToDataFrame('Metapopulations_example.csv')
graph = my_model.compositionPlot('Metapopulations_example.png', data, num_top_genomes=6, track_specific_genomes=['BADD'])
