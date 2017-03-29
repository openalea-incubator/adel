from openalea.core.alea import *
from openalea.mtg.traversal import pre_order2
from alinea.adel import mtg
import alinea.adel.fitting as fitting
from alinea.adel.symbol import build_symbols
from openalea.plantgl.all import *

g = None

def adel_mtg():
    global g
    if g is not None:
        return g
    pm = load_package_manager()
    g = run(('alinea.adel.tutorials','mtgparam'), {},pm=pm, vtx_id=5)

    g = g[0]
    

    return g

def leaves_db():
    import cPickle as Pickle
    fn = r'../adel/data/leaves_simple.db'
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves = fitting.fit_leaves(leaves, 9)
    
    db = leaves[0]
    functions = build_symbols(db)
    return functions

def test_adel_mtg():
    g = adel_mtg()
    #print g
    symbols = leaves_db()

    assert len(g) == len(g.sub_mtg(g.root))

    print symbols.keys()
    g=mtg.mtg_turtle(g,symbols)
    return g


def test_dynamic_mtg():
    """ Add to the each metamers a thermal time (e.g. 10).
    Length = length * (thermal_time - start_thermaltime)/(end-start)
    if the thermal time is lower that the start thermal time, 
    do not proceed more.
    """

    g = adel_mtg()

    v = g.component_roots_at_scale_iter(g.root, scale=3).next()
    tt = 0
    dtt = 10.
    for metamer in pre_order2(g, v):
        nm = g.node(metamer)
        for node in nm.components():
            node.start_tt = tt
            node.end_tt = tt+dtt
        tt += dtt

    return g, tt

def display():

    g, max_time = test_dynamic_mtg()
    symbols = leaves_db()

    for time in range(1,150):
        g.properties()['geometry'] = {}
        g = mtg.mtg_turtle_time(g, symbols,time)
        scene = g.to_plantgl()
        Viewer.display(scene)
        
    
