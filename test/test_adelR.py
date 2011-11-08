# Test AdelR
#
#
from alinea.adel.AdelR import devCsv,setAdel,RunAdel,genGeoLeaf,genGeoAxe
from alinea.adel.mtg import mtg_factory

def test_adelR(nplants,dd):
    devT = devCsv('./data/axeTCa0N.csv','./data/dimTCa0N.csv','./data/phenTCa0N.csv','./data/earTCa0N.csv','./data/ssi2sen.csv')
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT,geoLeaf,geoAxe,nplants)
    string = RunAdel(dd,pars)
    return pars,string

def test2():
    pars, d = test_adelR(1,500)
    g = mtg_factory(d)
    return g
