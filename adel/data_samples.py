""" Provides sample data for testing Adel without reference to file paths, nor needs to call package manager
"""

import os

datadir = os.path.dirname(__file__)

def leaves_db():
    import cPickle as Pickle
    import alinea.adel.fitting as fitting
    fn = datadir + '/data/leaves_simple.db'
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves,discard = fitting.fit_leaves(leaves, 9)
    return leaves
    
def devT():
    from alinea.adel.AdelR import devCsv
    
    axeT = datadir + '/data/axeTCa0N.csv'
    dimT = datadir + '/data/dimTCa0N.csv'
    phenT = datadir + '/data/phenTCa0N.csv'
    earT = datadir + '/data/earTCa0N.csv'
    ssisenT = datadir + '/data/ssi2sen.csv'
    
    return devCsv(axeT,dimT,phenT,earT,ssisenT)
    
def srdb():
    return datadir + '/data/SRSo.RData'
    
def xydb():
    return datadir + '/data/So99.RData'
    
def wheat_leaf_db():
    from alinea.adel.wheat.extract_wheat import extract_leaf_info
    import alinea.adel.fitting as fitting
    leaves= extract_leaf_info(xydb(),srdb())
    leaves,discard = fitting.fit_leaves(leaves, 9)
    return leaves