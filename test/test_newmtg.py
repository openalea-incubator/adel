from alinea.adel.newmtg import *
from alinea.adel.mtg_interpreter import *
from openalea.plantgl.all import *
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe
import alinea.adel.fitting as fitting


def adelR(nplants,dd):
    devT = devCsv('./data/axeTCa0N.csv','./data/dimTCa0N.csv','./data/phenTCa0N.csv','./data/earTCa0N.csv','./data/ssi2sen.csv')
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants)
    cantable = RunAdel(dd,pars)
    return pars,cantable
    
def leaves_db():
    import cPickle as Pickle
    fn = r'../adel/data/leaves_simple.db'
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves = fitting.fit_leaves(leaves, 9)
    return leaves
       
#from openalea.mtg import mtg2axialtree
#affichage : aller voir __str__dans mtg.py
#dp = {'plant':[1,1,1,1,1], 'axe' : [0,0,1,1,2], 'numphy': [1,2,1,2,1], 'Laz': [0,180,90,180,180],'Einc':[0,0,45,-45,45], 'El': [1,1,1,1,1],'Ev':[1,1,1,1,1], 'Esen':[0,0,0,0,0], 'Ed':[0.1,0.1,0.1,0.1,0.1]}
#dp = {'plant':[1,1,1,1,1], 'axe' : [0,0,1,1,2], 'numphy': [1,2,1,2,1], 'Laz': [0,180,90,180,0],'Einc':[0,0,45,-30,45], 'El': [1,1,1,1,1],'Ev':[1,1,1,1,1], 'Esen':[0,0,0,0,0], 'Ed':[0.1,0.1,0.1,0.1,0.1]}
#g=mtg_factory(dp,adel_metamer)
#g.display()

#more test
p,d= adelR(1,1000)
leaves=leaves_db()
g=mtg_factory(d,adel_metamer,leaf_db=leaves, stand=[((0,0,0),0),((10,0,0),90), ((0,10,0), 0)])
g=mtg_interpreter(g)
scene = plot3d(g)
Viewer.display(scene)
