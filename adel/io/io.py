# Standard import

import cPickle as Pickle

import numpy as np
from rpy2 import robjects
from rpy2.robjects.numpy2ri import numpy2ri

from alinea.adel.AdelR import RlistAsDict,readRData,saveRData,csvAsDict,dataframeAsdict,canL2canS
from alinea.adel.mtg import mtg_factory, duplicate, thermal_time
from openalea.mtg.io import lpy2mtg, mtg2lpy

def dataframe(d):
    """ convert a dict of numbers to an RDataframe  """
    df = {}
    if d is None:
        return robjects.r('as.null()')
    else:
        for k, v in d.iteritems():
            df[k] = numpy2ri(np.array(v))
    dataf = robjects.r['data.frame'](**df)
    return dataf



def load_leaf_data(fn):
    """  Load leaf data obtained by measurement. """ 
    leaves = {}
    try:
        f = open(fn)
        leaves = Pickle.load(f)
    finally:
        f.close()

    return leaves,

def to_canestra(can_scene):
    return can_scene.to_canestra()

def to_plantgl(can_scene, 
               leaf_material, 
               stem_material, 
               soil_material):
    return can_scene.to_plantgl(leaf_material, stem_material, soil_material)
