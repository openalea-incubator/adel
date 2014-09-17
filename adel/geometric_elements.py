"""
New proposal for computing organ shapes
"""

import pandas
import numpy

import alinea.adel.fitting as fitting
import alinea.adel.data_samples as adel_data
import alinea.adel.AdelR as adelR

# function for handling dynamical leaf shape databases    

class LeafGeometry(object):
    
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
        
    def leafdb(self):

        leaves = {}
        xy = self.xydb
        sr = self.srdb
        for k in xy.keys():
            leaves[k]=[]
            for i in range(len(xy[k])):
                leaves[k].append((xy[k][i][0], xy[k][i][1], sr[k][0], sr[k][1]))
                
        leaves,discard = fitting.fit_leaves(leaves, self.discretisation_level)
        
        return leaves

        
    
def dynamic_key(age, bins = [-10, 0.5, 1, 2, 3, 4, 10]):
    index = numpy.searchsorted(bins, age)
    if index ==0:
        index = 1 # age below first value are in firts interval
    return '(%s, %s]'%(str(bins[index-1]), str(bins[index]))
    
def leaf_trajectories(dfxy, dfsr, bins = [-10, 0.5, 1, 2, 3, 4, 10], ntraj = 10, tol_med = 0.1):
    """
    Return a dynamic leaf database compatible with adel and appropriate dynamic key selecting functions ({Lindex:{num_traj:{age_class:xy_dict}}})
    
    - dfxy and srdb are dataframe containing the data
    - bins is None if leaf are static or define age clas otherwise
    - ntraj and tol_med control  the number of trajectory to sample among leaves of a given Lindex and of given age that have mean_angle +/- tol*med
    
    """
                   
    import random
          
    dfxy['age'] = dfxy['HS'] - dfxy['rank'] + 1
    #cut intervalle de 0 a 1, etc.
    dfxy_cut = pandas.cut(dfxy.age, bins)
    dfxy['age_class'] = dfxy_cut
    
    # use mean_angle to filter / group leaves
    mean_angle = dfxy.groupby('inerv').apply(lambda x: numpy.mean(abs(numpy.arctan2(numpy.diff(x['y'].values),numpy.diff(x['x'].values)))))

    #filter leaves that are above/below med*tol
    def filter_leaves(x):
        angles = mean_angle[set(x['inerv'])]
        med = angles.median()
        valid_angles = angles[(angles >= (1 - tol_med) * med) & (angles <= (1 + tol_med) * med)]
        return x[x['inerv'].isin(set(valid_angles.index))]
    
    validxy = dfxy.groupby(('Lindex','age_class'), group_keys=False).apply(filter_leaves)
    grouped = validxy.groupby(('Lindex','age_class'))
    
    # build trajectories
    trajectories = {k:{} for k in set(validxy['Lindex'])}
    for i in range(ntraj):
        for k in set(validxy['Lindex']):
            trajectories[k][i] = {}
            for t in set(validxy['age_class']):
                x = grouped.get_group((k,t))
                trajectories[k][i][t] = x.ix[x['inerv'] == random.sample(set(x['inerv']),1),['x','y']].to_dict('list')
    
    srdb = {k:v.ix[:,['s','r']].to_dict('list') for k, v in dfsr.groupby('Lindex')}
    
    return trajectories, srdb
    
