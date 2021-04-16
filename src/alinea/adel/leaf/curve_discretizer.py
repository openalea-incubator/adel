#from openalea.plantgl import all as pgl

def curve_discretizer(curve, n=100):
    '''    
    '''
    x = []; y = []; 

    import numpy as np
    u = np.linspace(curve.firstKnot, curve.lastKnot, n)
    pts = list(map(curve.getPointAt, u))
    x,y = list(zip(*pts))

    # return outputs
    return x, y,

def points2curve(x,y):
    from openalea.plantgl import all as pgl
    points = list(zip(x,y))
    return pgl.NurbsCurve2D.fit(points)
