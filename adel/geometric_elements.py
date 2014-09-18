"""
New proposal for computing organ shapes
"""

import numpy

import alinea.adel.fitting as fitting
import alinea.adel.data_samples as adel_data
import alinea.adel.AdelR as adelR

# function for handling dynamical leaf shape databases    

class Leaves(object):
    
    def __init__(self, xydb = None, srdb = None, geoLeaf = None, dynamic = False, discretisation_level = 9):
    
        if xydb is None:
            xydb = adel_data.xydb()
        
        if srdb is None:
            srdb = adel_data.srdb()
            
        if geoLeaf is None:
            geoLeaf = adelR.genGeoLeaf()

        self.xydb = xydb
        self.srdb = srdb
        self.geoLeaf = geoLeaf
        self.dynamic = dynamic
        self.discretisation_level = discretisation_level
        
        leaves = {}
        xy = self.xydb
        sr = self.srdb
        for k in xy.keys():
            leaves[k]=[]
            for i in range(len(xy[k])):
                leaves[k].append((xy[k][i][0], xy[k][i][1], sr[k][0], sr[k][1]))        
        leaves,discard = fitting.fit_leaves(leaves, self.discretisation_level)
        self.leaves = leaves
        
    def get_leaf(self, lindex, lseed):
        key, index = adelR.leaf_keys(lindex, lseed, self.srdb)
        leaf = self.leaves[key][index]
        return leaf
    
def dynamic_key(age, bins = [-10, 0.5, 1, 2, 3, 4, 10]):
    index = numpy.searchsorted(bins, age)
    if index ==0:
        index = 1 # age below first value are in firts interval
    return '(%s, %s]'%(str(bins[index-1]), str(bins[index]))
    
