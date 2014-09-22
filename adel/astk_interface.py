""" Class instanciating a wheat canopy and complying to astk canopy interface
"""

try:
    import cPickle as pickle
except:
    import pickle
import os

import numpy


from alinea.adel.AdelR import setAdel,RunAdel,genGeoAxe, checkAxeDyn, getAxeT, getPhenT, getPhytoT
from alinea.adel.newmtg import *
import alinea.adel.data_samples as adel_data
from alinea.adel.mtg_interpreter import *
from alinea.astk.TimeControl import *
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand


def get_normal_dist(nb_plants=10, sigma=30.):
    """ Calculate the "best possible" distribution of emergence dates for a given number of plants 
    as a function of the standard deviation sigma in degree days. """
    N = 1000
    norm_list = numpy.random.normal(scale=sigma, size=N)
    h = numpy.histogram(norm_list, bins=nb_plants)
    classes = h[1]
    distri_plants = h[0]*nb_plants/float(N)
    round_distri = [numpy.round(d) for d in distri_plants]
    missing_pl = nb_plants - sum(round_distri)
    while missing_pl > 0:
        gap = round_distri - distri_plants
        round_distri[np.argmin(gap)] += 1
        missing_pl = nb_plants - sum(round_distri)
    while missing_pl < 0:
        gap = distri_plants - round_distri
        round_distri[np.argmin(gap)] -= 1
        missing_pl = nb_plants - sum(round_distri)
    return numpy.hstack([numpy.linspace(h[1][i], h[1][i+1], d+2)[1:-1] for i, d in enumerate(round_distri) if d>0])

    
class AdelWheat(object):
    
    def __init__(self, nplants = 1, nsect = 1, devT = None, sample = 'random', seed = None, leaves = None, stand = None, aspect='square', convUnit = 0.01, thermal_time_model = None, incT=60, dinT=5, dep = 7, run_adel_pars = {'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'endLeaf1': 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6, 'HSstart_inclination_tiller': 1, 'rate_inclination_tiller': 30}, split=False, face_up=False, aborting_tiller_reduction = 1.0, classic=False, leaf_db = None, positions = None):
    
        if devT is None: 
            devT = adel_data.devT()
            
        if leaf_db is not None:
            print('!!!!Warning!!!! leaf_db argument is deprecated, use adel.geometric_elements.Leaves class instead')
        if positions is not None:
            print('!!!!Warning!!!! positions argument is deprecated, use stand = adel.Stand class instead')
            
        if leaves is None:
            leaves = Leaves()
            
        if stand is None:
            stand = AgronomicStand(sowing_density = 250, plant_density=250, inter_row=0.15)
            positions
            
        if thermal_time_model is None:
            thermal_time_model = DegreeDayModel(Tbase=0)

        geoAxe = genGeoAxe(incT=incT,dinT=dinT,dep=dep)
        self.pars = setAdel(devT,leaves.geoLeaf,geoAxe,nplants, seed = seed, sample=sample, xydb = leaves.xydb, srdb=leaves.srdb)
        self.leaves = leaves
        self.nsect = nsect
        self.nplants, self.domain, self.positions = stand.stand(nplants, aspect=aspect)
        self.domain_area = abs(self.domain[1][0] - self.domain[0][0]) * abs(self.domain[1][1] - self.domain[0][1]) * convUnit * convUnit
        self.thermal_time = thermal_time_model
        self.run_adel_pars = run_adel_pars
        self.split = split
        self.face_up = face_up
        self.aborting_tiller_reduction = aborting_tiller_reduction
        self.classic = classic
        self.convUnit = convUnit
    
    def timing(self, delay, steps, weather, start_date):
        """ compute timing and time_control_sets for a simulation between start and stop. 

        Return 0 when there is no rain
        """ 
        timestep = delay
        start_date= weather.str_to_datetime(start_date)
        temp = [d['temperature_air'] for d in weather.split_weather(timestep, start_date, steps)]
       
        return (TimeControlSet(Tair = temp[i / int(delay)], dt = delay) if not i % delay  else TimeControlSet(dt=0) for i in range(steps))

    
    def setup_canopy(self, age = 10):
    
        self.canopy_age = age
        canopy = RunAdel(age, self.pars, adelpars=self.run_adel_pars)
        if self.positions is not None:
            stand = [(pos,0) for pos in self.positions]
        else:
            stand = None
        g = mtg_factory(canopy, adel_metamer, leaf_sectors=self.nsect, leaves=self.leaves, stand=stand, split=self.split, aborting_tiller_reduction=self.aborting_tiller_reduction)
        g = mtg_interpreter(g, self.leaves, face_up = self.face_up, classic= self.classic)
        return g

    def checkAxeDyn(self, dates, density=1):
        return checkAxeDyn(self.pars, dates, density)
        
    def axeT(self):
        return getAxeT(self.pars)
        
    def phenT(self, axe='MS'):
        return getPhenT(self.pars, axe=axe)
        
    def phytoT(self, axe='MS'):
        return getPhytoT(self.pars, axe=axe)
        
    def grow(self, g, time_control):
    
        try: #old interface
            if time_control.dt <= 0:
                dday=0.
            else:
                dday = time_control.Tair.mean()
        except:
            data = time_control
            tt = self.thermal_time(data.index, data)
            dday = tt[-1]
        
        #refg = self.setup_canopy(age = self.canopy_age)
        self.canopy_age += dday
        newg = self.setup_canopy(age = self.canopy_age)
        #newg = mtg_update(newg, g, refg)
        move_properties(g, newg)
        
        return newg
   
    def grow_dd(self, g, dday):
        #refg = self.setup_canopy(age = self.canopy_age)
        self.canopy_age += dday
        newg = self.setup_canopy(age = self.canopy_age)
        #newg = mtg_update(newg, g, refg)
        move_properties(g, newg)
        return newg
   
    def plot(self, g):
        from openalea.plantgl.all import Viewer
        s = plot3d(g)
        Viewer.display(s)
        return s
        
    def scene(self, g):
        return plot3d(g)
        
    def get_exposed_areas(self, g, convert=False, TT=None):
        areas = exposed_areas(g)
        if convert:
            areas = exposed_areas2canS(areas)
        if TT is None:
            TT = self.canopy_age
        areas['TT'] = TT
        return areas
        
    def save(self, g, index = 0, dir = './adel_saved'):
        if not os.path.exists(dir):
            os.mkdir(dir)
        s = self.scene(g)
        geom = {sh.id:sh.geometry for sh in s}
        g.remove_property('geometry')
        fgeom = dir + '/scene%04d.bgeom'%(index)
        fg = dir+'/adel%04d.pckl'%(index)
        s.save(fgeom, 'BGEOM')
        f = open(fg, 'w')
        pickle.dump([g,self.canopy_age], f)
        f.close()
        #restore geometry
        g.add_property('geometry')
        g.property('geometry').update(geom)
        return fgeom,fg
        
    def load(self, index=0, dir = './adel_saved'):
        from openalea.plantgl.all import Scene
        
        fgeom = dir + '/scene%04d.bgeom'%(index)
        fg = dir + '/adel%04d.pckl'%(index)
        
        s = Scene()
        s.read(fgeom, 'BGEOM')
        geom = {sh.id:sh.geometry for sh in s}
    
        f = open(fg)
        g, TT = pickle.load(f)
        f.close()
        
        self.canopy_age = TT
        
        g.add_property('geometry')
        g.property('geometry').update(geom)
        
        return g, TT
        

def adelwheat_node(nplants = 1, positions = None, nsect = 1, devT = None, leaf_db = None, sample = 'random', seed = None, thermal_time_model = None):
    model = AdelWheat(nplants = nplants, positions = positions, nsect = nsect, devT = devT, leaf_db = leaf_db, sample = sample, seed = seed, thermal_time_model = thermal_time_model)
    return model
    
        
#user friendly macros
from alinea.adel.stand.stand import agronomicplot
from alinea.astk.plant_interface import *

def initialise_stand(age=0., length=0.1, width=0.2, sowing_density=150, 
                     plant_density=150, inter_row=0.12, nsect=1, seed = None, sample='random'):
    """ Initialize a wheat canopy.
    
    Parameters
    ----------
    age: float
        Age of the canopy at initialization (in degree days)
    length: float
        Plot dimension along row direction (in m)
    width: float
        Plot dimension perpendicular to row direction (in m)
    sowing density: int
        Density of seeds sawn (in seeds.m-2)
    plant_density: int
        Density of plants that are present (after loss due to bad emergence, 
        early death...) (in plants.m-2)
    inter_row: float
        Distance between rows (in m)
    seed: float
        random seed used by adel to sample plants in the data
    sample : string
        type of sampling. 'random' or 'sequence'
    
    Returns
    -------
    g: MTG
        Wheat canopy
    wheat: instance
        Wheat instance of AdelWheat
    domain_area: float
        Soil surface occupied by plants (inverse of density) (in m2)
    domain : tuple
        tuple of coordinates defining the domain

    """
    nplants, positions, domain, domain_area, convUnit = agronomicplot(length=length, 
                                                            width=width, 
                                                            sowing_density=sowing_density, 
                                                            plant_density=plant_density,
                                                            inter_row=inter_row)
    wheat = AdelWheat(nplants=nplants, positions = positions, nsect=nsect, seed= seed, sample=sample)
    g,_ = new_canopy(wheat,age=age)
    return g, wheat, domain_area, domain, convUnit
        