"""
New proposal for computing organ shapes
"""

import numpy
from scipy.integrate import simps

import openalea.plantgl.all as pgl
from math import degrees, radians, pi, cos, sin

from alinea.adel.exception import *
import alinea.adel.fitting as fitting
import alinea.adel.data_samples as adel_data
import alinea.adel.AdelR as adelR

# function for handling dynamical leaf shape databases    
def incline_leaf(shape, inclin, relative_angle = True):
    """ transform a xysr tuple representing leaf shape to get a given angle at leaf base.
     - angle the desired angle (deg)
     - if relative_angle == True, angle is interpreted as a multiplier to original shape angle
     """   
    Linc = inclin
    x, y = shape[0], shape[1]
    init_angle = pgl.angle((x[1]-x[0], y[1]-y[0]),(0,1))

    if relative_angle:
        angle = Linc * init_angle
        angle = min(pi, angle)
    else:
        angle = radians(Linc)
    
    rotation_angle = init_angle - angle

    # rotation of the midrib
    cos_a = cos(rotation_angle); sin_a = sin(rotation_angle)

    x1 = x[0] + cos_a*x - sin_a*y
    y1 = y[0] + sin_a*x + cos_a*y
    leaf= x1, y1, shape[2], shape[3]
    
    return leaf


class Leaves(object):
    
    def __init__(self, xydb = None, srdb = None, geoLeaf = None, dynamic_bins = None, discretisation_level = 9, twist = 0):
    
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
        self.twist = twist
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
        if isinstance(leaf, dict):
            leaf = leaf['x'], leaf['y'], leaf['s'], leaf['r']
        return leaf
        
    def get_sr(self, leaf_key):
        key, index, age_index = leaf_key
        sr = self.srdb[key]
        if isinstance(sr, dict):
            sr = sr['s'], sr['r']
        return sr
        
    def blade_elt_area(self, leaf_key, Lshape, Lwshape, sr_base, sr_top):
        """ surface of a blade element, positioned with two relative curvilinear absisca"""
       
        S=0
        sr_base = min([1,max([0,sr_base])])
        sr_top = min([1,max([sr_base,sr_top])])
        if leaf_key is not None:
            s,r = self.get_sr(leaf_key)
            sre = [sr for sr in zip(s,r) if (sr_base < sr[0] < sr_top)]
            if len(sre) > 0:
                se,re = zip(*sre)
                snew = [sr_base] + list(se) + [sr_top]
                rnew = [numpy.interp(sr_base,s,r)] + list(re) + [numpy.interp(sr_top,s,r)]
            else:
                snew = [sr_base, sr_top]
                rnew = [numpy.interp(sr_base,s,r), numpy.interp(sr_top,s,r)]

            S = simps(rnew,snew) * Lshape * Lwshape

        return S

    def mesh(self, leaf_key, L_shape, Lw_shape, length, s_base, s_top, incline=1, flipx = False):
        """ Compute mesh for a leaf element.
            - shape is a x,y,s,r tuple descriibing leaf shape
            - L_shape is the length of the scaled shape
            - Lw_shape is the width of the scaled shape
            - length is the total visible length to be meshed
            - s_base and s_top are relative proportion (on length) of the element to represent
        """
        shape = self.get_leaf(leaf_key)

        shape = incline_leaf(shape, incline)
        if flipx:
            shape = (-shape[0],)+shape[1:] # to position leaves along tiller emitted
        leaf_mesh = fitting.mesh4(shape, L_shape, length, s_base, s_top, Lw_shape, twist=self.twist)
        if leaf_mesh:
            pts, ind = leaf_mesh
            if len(ind) < 1:
                raise  AdelError('ERROR less than 1 triangles')
                mesh = None
            else:
                mesh = fitting.plantgl_shape(pts, ind)
        else:
            if length > 0:
                print 'ERROR No mesh', s_base, s_top, length
                pass
            mesh = None

        return mesh
    
    def form_factor(self):
        """
        return form factor for each key in sr_db
        """
        
        return {k:simps(self.srdb[k]['r'], self.srdb[k]['s']) for k in self.srdb}
        

