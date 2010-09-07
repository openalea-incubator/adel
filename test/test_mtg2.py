from alinea.adel.mtg import *

def test1():
    fn = r'data/canopydesc.csv'
    g = topological_table_to_mtg(fn)
    assert g.nb_scales() == 5
    assert g.max_scale()==4
    assert len(g) == 93
    
    return g
