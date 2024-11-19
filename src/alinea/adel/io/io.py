# Standard import

import pickle as Pickle

import numpy
from rpy2 import robjects

r = robjects.r

# needed for wralea
from alinea.adel.AdelR import (
    RlistAsDict,
    readRData,
    saveRData,
    csvAsDict,
    dataframeAsdict,
    canL2canS,
    dataframe,
)
from alinea.adel.mtg import (
    mtg_factory,
    duplicate,
    thermal_time,
    apply_property,
    to_plantgl,
    to_canestra,
)
from openalea.mtg.io import lpy2mtg, mtg2lpy

import openalea.core.path as path


def load_leaf_data(fn):
    """Load leaf data obtained by measurement."""
    leaves = {}
    try:
        f = open(fn)
        leaves = Pickle.load(f)
    finally:
        f.close()

    return (leaves,)
