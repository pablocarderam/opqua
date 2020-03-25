#TODO: write example package usage

'''

'''

from opqua.model import Model

m = Model()
m.newSetup('s')
m.newPopulation('p1','s')

#TODO: test fitness functions
#TODO: test pathogens
#TODO: test interventions

results = m.run(0,10,"~/Desktop/test.csv")
