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
    
    def get_leaf_key(self, lindex, lseed, age=None):
    # to do return one default leaf even if key error occur ?
        key, index = adelR.leaf_keys(lindex, lseed, self.srdb)
        age_index = age
        if age is not None:
            age_index = numpy.searchsorted(self.bins, age)
            if age_index ==0:
                age_index = 1 # age below first value are in firts interval
            if age_index >= len(self.bins):
                age_index = len(self.bins) - 1# age above last value are set in last interval
            age_index = '(%s, %s]'%(str(self.bins[age_index-1]), str(self.bins[age_index]))
        return key, index, age_index
    
    def get_leaf(self, leaf_key):
        key, index, age_index = leaf_key
        if age_index is None:
            leaf = self.leaves[key][index]
        else:
            leaf = self.leaves[key][index][age_index]
        return leaf
        
    def blade_elt_area(self, leaf_key, Lshape, Lwshape, sr_base, sr_top):
        """ surface of a blade element, positioned with two relative curvilinear absisca"""
        
        from scipy.integrate import simps
        from numpy import interp, arange

        S=0
        sr_base = min([1,max([0,sr_base])])
        sr_top = min([1,max([sr_base,sr_top])])
        if leaf_key is not None:
            x,y,s,r = self.get_leaf(leaf_key)
            sre = [sr for sr in zip(s,r) if (sr_base < sr[0] < sr_top)]
            if len(sre) > 0:
                se,re = zip(*sre)
                snew = [sr_base] + list(se) + [sr_top]
                rnew = [interp(sr_base,s,r)] + list(re) + [interp(sr_top,s,r)]
            else:
                snew = [sr_base, sr_top]
                rnew = [interp(sr_base,s,r), interp(sr_top,s,r)]

            S = simps(rnew,snew) * Lshape * Lwshape

        return S

