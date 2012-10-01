from alinea.adel.newmtg import *
from openalea.mtg.io import display_tree

def display(g,scale=-1):
    if scale < 0:
        print g
    else:
        print '\n'.join(display_tree(g,scale))
       
#from openalea.mtg import mtg2axialtree
#affichage : aller voir __str__dans mtg.py
dp = {'plant':[1,1,1,1,2], 'axe' : [0,0,1,1,0], 'numphy': [1,2,1,2,1]}

g=mtg_factory(dp,adel_metamer)
print g