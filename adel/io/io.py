# Standard import

import cPickle as Pickle

from alinea.adel.AdelR import RlistAsDict,readRData,saveRData,csvAsDict

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
