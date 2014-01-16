from alinea.adel.astk_interface import AdelWheat
from alinea.astk.TimeControl import TimeControlSet
from alinea.adel import postprocessing

nplants=1
wheat = AdelWheat(nplants=nplants)

def test_usage():
    g = wheat.setup_canopy()
    timing = [TimeControlSet(dt=100) for i in range(2)]
    for tc in timing:
        wheat.plot(g)
        wheat.grow(g,tc)
    return g
    
def test_postprocessing(age=100):
    g = wheat.setup_canopy(age=age)
    out = wheat.get_exposed_areas(g,convert=True)
    domain_area = 1
    axis_df = postprocessing.axis_statistics(out, domain_area)
    plot_df = postprocessing.plot_statistics(axis_df, nplants, domain_area)
    return axis_df, plot_df 

