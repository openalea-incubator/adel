"""
Define generic symbols for the different organs.

Symbols are:
    LeafElement( tissue_type (int), 
                 final_length, length, radius_max, s_base, s_top ([0,1]),
                 leaf_id (rank), seed (random) )
    StemElement( tissue type (int) , length, diam_base, diam_top)
"""
# Define symbol classes used by the Turtle like
# Leaf, Stem, PStem, Apex, ApexR

import random, math
import openalea.plantgl.all as pgl
import fitting
import label
import cPickle as Pickle
from copy import deepcopy
from itertools import chain
from math import degrees, radians, pi, cos, sin
import numpy

classic = False

def optical(tissue_type):
    if tissue_type in (1,3,5,7,9):
        opt = 1
    elif tissue_type in (2,4,6,8,10):
        opt = 2
    return int(opt)


class Symbol(object):
    def _mesh(self):
        pass

    def _can_label(self, leaf_number, optical_species ):
        _label = label.Label()
        _label.leaf_id = leaf_number
        _label.optical_id = optical_species
        return _label


# Christian: add one arg( length ) in the __call.
# functions are commented.
# s_base and s_top are relative to current length and not final_length.

class LeafElement(Symbol):
    def __init__(self, database, seed = None, relative_angle=True):
        self.database = database
        self.seed = seed
        self.relative_angle = relative_angle

    def __call__( self, tissue_type, final_length, length, radius_max, s_base, s_top, leaf_rank, seed=1, *args ):
    #def __call__( self, optical_species, final_length, radius_max, s_base, s_top, leaf_rank, seed ):

        leaf_rank = max(1, leaf_rank% 999)
        opt = optical(tissue_type)
        res = {}

        
        res['geometry'] = self._mesh(leaf_rank, seed, final_length, length, s_base, s_top, radius_max, *args)
        #res['geometry'] = self._mesh(leaf_rank, seed, final_length, final_length, s_base, s_top, radius_max)
        res['label'] = self._can_label(leaf_rank, opt)
        res['tissue_type'] = tissue_type
        return res


    def _mesh(self, leaf_rank, seed, total_length, length, s_base, s_top, radius_max, *args):
        
        db = self.database
        rank_max = max(map(int,db.keys()))
        rank = leaf_rank
        rank = min(rank, rank_max)
        #choisi la liste de leaves du rang, ou rang + 1 si clef absente ourag -1 ou liste vide sinon
        leaves = db.get(str(rank), db.get(str(rank+1), db.get(str(rank-1), [])))
        n = len(leaves)
        if n == 0:
            raise "Not enough leaves data at rank %d"%rank

        if self.seed is None:
            random.seed(seed)
        else: 
            random.seed(self.seed)

        if args and len(args) > 1 and args[1] >= 0:
            i = args[1]
        else:
            i = random.randint(0,n-1)
            
        leaf = leaves[i]

        # Rotation of the midrib of the leaf to set the insertion angle a relative fraction of the angle (reference beeing the  vertical)
        if args and args[0] >= 0:
            Linc = args[0]
            x, y = leaf[0], leaf[1]
            init_angle = pgl.angle((x[1]-x[0], y[1]-y[0]),(0,1))

            if self.relative_angle:
                angle = Linc * init_angle
                angle = min(math.pi, angle)
            else:
                angle = math.radians(Linc)
            
            rotation_angle = init_angle - angle

            # rotation of the midrib
            cos_a = cos(rotation_angle); sin_a = sin(rotation_angle)

            x1 = x[0] + cos_a*x - sin_a*y
            y1 = y[0] + sin_a*x + cos_a*y

            leaf = (x1, y1) + leaf[2:]

        leaf_mesh = fitting.mesh4(leaf, total_length, length, s_base, s_top, radius_max)
        if leaf_mesh:
            pts, ind = leaf_mesh
            if len(ind) < 2:
                mesh = None
            else:
                mesh = fitting.plantgl_shape(pts, ind)
        else:
            mesh = None

        return mesh


class StemElement(Symbol):
    def __init__(self, database, seed=None):
        #self.database = database
        self.database = []
        self.tessel = pgl.Tesselator()
        self.seed = seed

    def __call__( self, tissue_type, length, diameter_base, diameter_top, *args ):
        leaf_rank = 0
        opt = optical(tissue_type)
        res = {}
        res['geometry'] = self._mesh(length, diameter_base, diameter_top)
        res['label'] = self._can_label(leaf_rank, opt)
        res['tissue_type'] = tissue_type
        return res

    def _mesh(self, length, diameter_base, diameter_top):
        if self.seed is not None or classic:
            solid = True
            slices = 3
            diameter = diameter_base
            stem = pgl.Tapered(diameter_base/2., diameter_top/2., pgl.Cylinder(1., length , solid, slices)) 
            stem.apply(self.tessel)
            mesh = self.tessel.triangulation
        else:
            mesh = slim_cylinder(length, diameter_base /2., diameter_top /2.)

        return mesh


def build_symbols(leaves, seed=None, relative_angle=True):

    symbols = {}
    symbols['LeafElement'] = LeafElement(leaves,seed=seed, relative_angle=relative_angle)
    symbols['StemElement'] = StemElement(leaves,seed=seed)
    return symbols


# geometric functions
def slim_cylinder(length, radius_base, radius_top):
    " Try to construct a cylinder with a low number of triangles. "
    pi = math.pi
    cos = math.cos
    sin = math.sin
    rb, rt = radius_base, radius_top
    a1, a2, a3 = 0, 2*pi/3., 4*pi/3.
    r = rb
    p1 = (r*cos(a1), r*sin(a1),0)
    p2 = (r*cos(a2), r*sin(a2),0)
    p3 = (r*cos(a3), r*sin(a3),0)
    r = rt
    q1 = (r*cos(a1+pi), r*sin(a1+pi),length)
    q2 = (r*cos(a2+pi), r*sin(a2+pi),length)
    q3 = (r*cos(a3+pi), r*sin(a3+pi),length)
    set = pgl.TriangleSet([p1, p2, p3, q1, q2, q3],
                      [(2,1,0), (3,4,5), (0,5,4), (0,4,2), (2,4,3), (3,1,2), (1,3,5), (5,0,1)])
    return set


def normal( set ):
    pts = set.pointList
    indices = set.indexList
    geoms = []
    for index in indices:
        i1, i2, i3 = index
        p1, p2, p3 = pts[i1], pts[i2], pts[i3]
        pt_bary = (p1+p2+p3)/3.
        n = (p2-p1) ^ (p3-p1)
        geoms.append(TriangleSet([p1, p2, p3], [(0,1,2)]))
        geoms.append(Polyline([pt_bary, pt_bary+ n]))
    return geoms

