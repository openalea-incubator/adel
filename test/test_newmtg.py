from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe


def adelR(nplants,dd):
    devT = devCsv('./data/axeTCa0N.csv','./data/dimTCa0N.csv','./data/phenTCa0N.csv','./data/earTCa0N.csv','./data/ssi2sen.csv')
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants)
    cantable = RunAdel(dd,pars)
    return pars,cantable
    

       
#from openalea.mtg import mtg2axialtree
#affichage : aller voir __str__dans mtg.py
dp = {'plant':[1,1,1,1,2], 'axe' : [0,0,1,1,0], 'numphy': [1,2,1,2,1]}

g=mtg_factory(dp,adel_metamer)
#g.display()

#more test
p,d= adelR(1,1000)
g=mtg_factory(d,adel_metamer)
g=mtg_interpreter(g)
scene = plot3d(g)
Viewer.display(scene)
