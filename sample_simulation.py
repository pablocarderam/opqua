#TODO: write example package usage

'''

'''

from opqua.model import Model

m = Model()
m.newSetup('s')
m.newPopulation('p1','s')
m.addPathogens( 'p1',{'AAAAAAAAAA':1} )
#TODO: test fitness functions
#TODO: test interventions

res = m.run(0,10,"~/Desktop/test.csv")
