from openalea.core.alea import *
from alinea.adel import mtg
import alinea.adel.fitting as fitting
from alinea.adel.symbol import build_symbols

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
    
    db = leaves
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
