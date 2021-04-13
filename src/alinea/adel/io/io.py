# Standard import

import pickle as Pickle

import numpy
from rpy2 import robjects
r = robjects.r

# needed for wralea
from alinea.adel.AdelR import RlistAsDict,readRData,saveRData,csvAsDict,dataframeAsdict,canL2canS
from alinea.adel.mtg import (mtg_factory, duplicate, thermal_time,
                             apply_property, to_plantgl, to_canestra)
from openalea.mtg.io import lpy2mtg, mtg2lpy

import openalea.core.path as path

def _is_iterable(x):
    try:
        x = iter(x)
    except TypeError:
        return False
    return True


def dataframe(d):
    """ convert a dict of numbers to an RDataframe  """
    df = {}
    if d is None:
        return r('as.null()')
    else:
        for k, v in d.items():
            rval = numpy.array(v)
            if not _is_iterable(v):
                v = [v]
            if 'NA' in numpy.array(v).tolist():
                df[k] = r['as.numeric'](rval)
            else :
                df[k] = rval
    dataf = r['data.frame'](**df)
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


