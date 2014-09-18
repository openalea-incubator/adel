"""
New proposal for computing organ shapes
"""

import numpy

import alinea.adel.fitting as fitting
import alinea.adel.data_samples as adel_data
import alinea.adel.AdelR as adelR

# function for handling dynamical leaf shape databases    
    

class Leaves(object):
    
    def __init__(self, xydb = None, srdb = None, geoLeaf = None, dynamic_bins = None, discretisation_level = 9):
    
        if xydb is None:
            xydb = adel_data.xydb()
        
        if srdb is None:
            srdb = adel_data.srdb()
            
        if geoLeaf is None:
            geoLeaf = adelR.genGeoLeaf()
                
        if dynamic_bins is None:
            self.dynamic = False
        else:
            self.dynamic = True

        self.xydb = xydb
        self.srdb = srdb
        self.geoLeaf = geoLeaf
        self.bins = dynamic_bins
        self.discretisation_level = discretisation_level
        self.fit_leaves()

    def fit_leaves(self):
        leaves = {}
        xy = self.xydb
        sr = self.srdb
        for k in xy.keys():
            leaves[k]=[]
            for i in range(len(xy[k])):
                if self.dynamic :
                    xysr = {age:(xy[k][i][age]['x'], xy[k][i][age]['y'], sr[k]['s'], sr[k]['r']) for age in xy[k][i]}
                else:
                    xysr = (xy[k][i][0], xy[k][i][1], sr[k][0], sr[k][1])
                leaves[k].append(xysr)
        leaves,discard = fitting.fit_leaves(leaves, self.discretisation_level, self.dynamic)
            
        self.leaves = leaves
    
    def get_leaf(self, lindex, lseed):
        key, index = adelR.leaf_keys(lindex, lseed, self.srdb)
        leaf = self.leaves[key][index]
        return leaf
    
    def dynamic_key(self, age):
        index = numpy.searchsorted(self.bins, age)
        if index ==0:
            index = 1 # age below first value are in firts interval
        return '(%s, %s]'%(str(self.bins[index-1]), str(self.bins[index]))
    
