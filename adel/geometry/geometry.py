from alinea.adel import fitting, symbol
from alinea.adel.mtg import mtg_turtle, mtg_turtle_time

def leaf_to_mesh(leaf, lmax, l, rmax):
    pts, ind = fitting.mesh3(leaf, lmax, l, rmax)
    mesh = fitting.plantgl_shape(pts, ind)
    return mesh

def leaf_element(leaf, lmax, l, s_base, s_top,  rmax):
    pts, ind = fitting.mesh4(leaf, lmax, l, s_base, s_top, rmax)
    mesh = fitting.plantgl_shape(pts, ind)
    return mesh

def symbols( database, seed=None ):
    return symbol.build_symbols(database, seed)

