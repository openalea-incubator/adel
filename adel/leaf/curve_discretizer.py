import numpy as np
from openalea.plantgl import all as pgl

def curve_discretizer(curve, n=100):
    '''    
    '''
    x = []; y = []; 

    u = np.linspace(curve.firstKnot, curve.lastKnot, n)
    pts = map(curve.getPointAt, u)
    x,y = zip(*pts)

    # return outputs
    return x, y,

def points2curve(x,y):
    points = zip(x,y)
    return pgl.NurbsCurve2D.fit(points)
