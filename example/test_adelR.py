# Test AdelR
#
#
from openalea.plantgl.all import Viewer
from pandas import DataFrame
import alinea.adel.data_samples as adel_data

leafdb = adel_data.wheat_leaf_db()

def adelR(nplants,dd):
    from alinea.adel.AdelR import setAdel,RunAdel,genGeoLeaf,genGeoAxe
    devT = adel_data.devT()
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants,xydb=adel_data.xydb(),srdb=adel_data.srdb())
    cantable = RunAdel(dd,pars)
    return pars,cantable
   
def getString(d):
    from alinea.adel.AdelR import dataframe, genString
    return genString(dataframe(d))
    
def canMTG(string):
    from alinea.adel.symbol import build_symbols
    from alinea.adel.mtg import CanMTG
    symbols = build_symbols(leafdb)
    g = CanMTG(symbols, string)
    s = g.to_plantgl()
    return g,s[0]
    
def newmtg(canopy,dec=10,az=0):
    from alinea.adel.newmtg import mtg_factory, adel_metamer
    from alinea.adel.mtg_interpreter import mtg_interpreter, plot3d
    g=mtg_factory(canopy,adel_metamer, leaf_sectors=1,leaf_db=leafdb, stand = [((dec,0,0),az)])
    g=mtg_interpreter(g)
    s = plot3d(g)
    return g,s

def test_positioning():
    d = {'plant':[1,1,1],'axe_id':['MS','T1','T1'],'ms_insertion':[0,1,1],'numphy':[1,1,2],
         'Laz': [0,90,180], 'Ll' :[3,3,3], 'Lv' :[3,3,3] , 'Lsen':[0,0,1], 'L_shape':[3,3,3], 'Lw_shape':[.3,.3,.3], 'Linc':[0,0,0],
         'Einc':[0,45,-45],'El':[1,1,1],'Ev':[1,1,1],'Esen':[0,0,0],'Ed': [0.1,0.1,0.1], 
         'Ginc':[0,0,0],'Gl':[0,0,0],'Gv':[0,0,0],'Gsen':[0,0,0],'Gd': [0.1,0.1,0.1], 
         'Epo':[1,1,1], 'Epos':[1,1,1], 'Gpo':[1,1,1], 'Gpos':[1,1,1], 'Lpo':[1,1,1], 'Lpos':[1,1,1],
         'LsenShrink':[1,1,1], 'LcType':[1,1,1], 'LcIndex':[1,1,1], 'Lr':[0,0,0]}
    chn = getString(d)
    g1,s1 = canMTG(chn)
    g2,s2 = newmtg(d,dec=5)
    Viewer.display(s1+s2)

     
   
def test(date=500,dec=10):
    pars, d = adelR(1,date)
    chn = getString(d)
    g1,s1 = canMTG(chn)
    g2,s2 = newmtg(d,dec=dec)
    Viewer.display(s1+s2)
    
  
#df=DataFrame(d)
    #df.Lv
    #df.ix[1:7,6:9]
    #df[['Gl','Gv','Ll','Lv','Lr','L_shape']]
# import numpy as np
# from alinea.popdrops.Rain import get_area_and_normal
# a1,_=get_area_and_normal(g1.property('geometry'))
# a1 = dict(((k,np.sum(v)) for k,v in a1.iteritems()))