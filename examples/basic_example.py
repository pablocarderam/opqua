
'''
Simple model
'''

from opqua.model import Model

my_model = Model()
my_model.newSetup('my_setup',default='vector-borne') # uses default parameters
my_model.newPopulation('my_population','my_setup')
my_model.addPathogensToHosts( 'my_population',{'AAAAAAAAAA':20} )
my_model.run(0,100)
data = my_model.saveToDataFrame('Basic_example.csv')
graph = my_model.compositionPlot('Basic_example.png', data, num_top_genomes=6)
