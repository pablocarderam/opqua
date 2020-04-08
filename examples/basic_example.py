
'''
Simple vector-borne model with susceptible and infected hosts and vectors.
'''

from opqua.model import Model

my_model = Model()
my_model.newSetup('my_setup',preset='vector-borne')
    # uses default parameters for vector-borne model
my_model.newPopulation('my_population','my_setup')
my_model.addPathogensToHosts( 'my_population',{'AAAAAAAAAA':20} )
my_model.run(0,200)
data = my_model.saveToDataFrame('Basic_example.csv')
graph = my_model.compartmentPlot('Basic_example.png', data)
