from alinea.adel import mtg

def CanMTGPlanter(g, positions, random_seed=0, azimuths = None):
    '''    arrange plant in a stand
    '''
    if hasattr(g, 'planter'):
        g.planter(positions, random_seed,azimuths)
    else:
        planter(g, positions, random_seed)
    return g,

