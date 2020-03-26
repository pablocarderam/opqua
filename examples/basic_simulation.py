#TODO: write example package usage

'''

'''

from opqua.model import Model

m = Model()
m.newSetup('s')
m.newPopulation('p1','s')
m.addPathogensToHosts( 'p1',{'AAAAAAAAAA':4} )
print(m.populations['p1'].num_infected_hosts)
#TODO: test fitness functions
#TODO: test interventions

res = m.run(0,10,"~/Desktop/test.csv")
