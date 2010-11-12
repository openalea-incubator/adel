from alinea.adel import mtg

def CanMTGPlanter(g, positions, random_seed=0):
    '''    arrange plant in a stand
    '''
    if hasattr(g, 'planter'):
        g.planter(positions, random_seed)
    else:
        planter(g, positions, random_seed)
    return g,

