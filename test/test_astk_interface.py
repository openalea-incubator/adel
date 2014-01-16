from alinea.adel.astk_interface import AdelWheat
from alinea.astk.TimeControl import TimeControlSet

wheat = AdelWheat()

def test_usage():
    g = wheat.setup_canopy()
    timing = [TimeControlSet(dt=100) for i in range(2)]
    for tc in timing:
        wheat.plot(g)
        wheat.grow(g,tc)
    return g
    
def test_postprocessing(age=10, convert=False):
    g = wheat.setup_canopy(age=age)
    out = wheat.get_exposed_areas(g,convert=convert)
    return g, out