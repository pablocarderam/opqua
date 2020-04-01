
'''
Simple model
'''

from opqua.model import Model

my_model = Model()
my_model.newSetup('my_setup',default="vector-borne") # uses default parameters
my_model.newPopulation('my_population','my_setup')
my_model.addPathogensToHosts( 'my_population',{'AAAAAAAAAA':20} )
my_model.run(0,100)
dat = my_model.saveToDataFrame("basic_example.csv")
