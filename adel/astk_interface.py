""" Class instanciating a wheat canopy and complying to astk canopy interface
"""

try:
    import cPickle as pickle
except:
    import pickle
import os

import numpy


from alinea.adel.AdelR import setAdel,RunAdel,genGeoAxe, checkAxeDyn, getAxeT, getPhenT, getPhytoT, saveRData
from alinea.adel.newmtg import *
import alinea.adel.data_samples as adel_data
from alinea.adel.mtg_interpreter import *
from alinea.astk.TimeControl import *
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.postprocessing import axis_statistics, plot_statistics
    
from openalea.mtg.algo import union
import openalea.plantgl.all as pgl

def replicate(g, target=1):

    current = g.nb_vertices(scale=1)
    if target <= current:
        return g
    
    g0 = g.sub_mtg(g.root)
    n0 = g.nb_vertices(scale=1)
    
    
    
    nduplication = int(numpy.log2(1. * target / current))   
    missing = target - current * numpy.power(2,nduplication)
    if missing >= n0:
        g_add = g.sub_mtg(g.root)
        missing -= n0
        
    for i in range(nduplication):
        g1 = g.sub_mtg(g.root)
        g = union(g, g1)
        n = g.nb_vertices(scale=1)
        if missing >= n and (missing - n) % n0 == 0:
            g1 = g.sub_mtg(g.root)
            g_add = union(g_add,g1)
            missing -= n
            
    g = union(g,g_add)
    
    if missing > 0:
        assert missing % n0 == 0
        nadd = int(missing / n0)
        for i in range(nadd):    
            g1 = g0.sub_mtg(g0.root)
            g = union(g,g1)

    return g

def transform_geom(geom, translation, rotation):
    if isinstance(geom, pgl.Geometry):
        geom = pgl.Translated(translation, pgl.AxisRotated((0,0,1),rotation, geom))
    elif isinstance(geom, pgl.Shape):
        geom = pgl.Shape(pgl.Translated(translation, pgl.AxisRotated((0,0,1),rotation, geom.geometry)))
    return geom

    
class AdelWheat(object):
    
    def __init__(self, nplants = 1, duplicate=1, nsect = 1, devT = None, sample = 'random', seed = None, leaves = None, stand = None, aspect='square', convUnit = 0.01, thermal_time_model = None, incT=60, dinT=5, dep = 7, run_adel_pars = {'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'endLeaf1': 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6, 'HSstart_inclination_tiller': 1, 'rate_inclination_tiller': 30, 'drop_empty':True}, split=False, face_up=False, aborting_tiller_reduction = 1.0, classic=False, leaf_db = None, positions = None, ssipars=None):
    
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
        self.stand=stand
            
        if thermal_time_model is None:
            thermal_time_model = DegreeDayModel(Tbase=0)

        geoAxe = genGeoAxe(incT=incT,dinT=dinT,dep=dep)
        
        self.devT = devT
        
        if self.stand.density_curve is None:
            self.nplants, self.domain, self.positions, area_m2 = stand.stand(nplants * duplicate, aspect=aspect)
        else:
            self.nplants, self.domain, self.positions, area_m2 = stand.smart_stand(nplants * duplicate)

        self.nplants_base = max(1,(self.nplants - self.nplants % duplicate) / duplicate)
        if self.nplants % self.nplants_base > 0:
            self.nplants_base = 1
        self.domain_area = area_m2
        self.plant_azimuths = numpy.random.random(self.nplants) * 360
        
        self.pars = setAdel(devT,leaves.geoLeaf,geoAxe,self.nplants_base, seed = seed, sample=sample, xydb = leaves.xydb, srdb=leaves.srdb, ssipars=ssipars)
        self.leaves = leaves
        self.nsect = nsect

        
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
        if self.stand.density_curve is not None:
            self.nplants, self.domain, self.positions, self.domain_area = self.stand.smart_stand(self.nplants, at=age)

        if self.positions is not None and self.nplants <= self.nplants_base:
            stand = zip(self.positions,self.plant_azimuths)
        else:
            stand = None
        g = mtg_factory(canopy, adel_metamer, leaf_sectors=self.nsect, leaves=self.leaves, stand=stand, split=self.split, aborting_tiller_reduction=self.aborting_tiller_reduction)
        g = mtg_interpreter(g, self.leaves, face_up = self.face_up, classic= self.classic)
        if self.nplants > self.nplants_base:
            g = replicate(g, self.nplants)
            pos = g.property('position')
            az = g.property('azimuth')
            lab = g.property('label')
            geom = g.property('geometry')
            for i,vid in enumerate(g.vertices(1)):
                lab[vid] = 'plant' + str(i+1)
                pos[vid] = self.positions[i]
                az[vid] = self.plant_azimuths[i]
                for gid in g.components_at_scale(vid, g.max_scale()):
                    if gid in geom:
                        geom[gid] = transform_geom(geom[gid], self.positions[i], self.plant_azimuths[i] * 1. / 180 *numpy.pi)
        
        return g

    def checkAxeDyn(self, dates=range(0,2000,100), density=None):
        if density is None:
            density=self.stand.plant_density
        return checkAxeDyn(self.pars, dates, density)
     
    def check_nff(self):
        """ Count the probability of occurence of MS with given nff
        """
        ms = numpy.array(self.devT['axeT']['id_axis']) == 'MS'
        nffs = numpy.array(self.devT['axeT']['N_phytomer'])[ms]
        counts = {n:nffs.tolist().count(n) for n in set(nffs)}
        probas = {n:nffs.tolist().count(n) * 1.0 / len(nffs) for n in set(nffs)}
        return counts, probas
        
        
    def check_primary_tillers(self):
        """ Count/estimate probabilitie of occurence of primary tillers
        """
        import re 
        axis = self.devT['axeT']['id_axis']
        ms = [e for e in axis if re.match('MS',e)]
        tillers = [e for e in axis if re.match('T.$',e)]
        counts = {n:tillers.count(n) for n in set(tillers)}
        probas = {n:tillers.count(n) * 1.0 / len(ms) for n in set(tillers)}
        return counts, probas
         
     
    def axeT(self):
        return getAxeT(self.pars)
        
    def phenT(self, axe='MS'):
        return getPhenT(self.pars, axe=axe)
        
    def phytoT(self, axe='MS'):
        return getPhytoT(self.pars, axe=axe)
        
    def save_pars(self):
        saveRData(self.pars,'plants','adel_pars.RData')
        
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
       
    def get_axis(self, g, plant='plant1', axe = 'MS'):
        """ return a new mtg representing an axe
        """
        p = [vid for vid in g.vertices(scale=1) if g.label(vid) == plant][0]
        ax = [vid for vid in g.components(p) if g.label(vid) == axe][0]
        return g.sub_mtg(ax,copy=True)
   
    def plot(self, g, property=None):
        from openalea.plantgl.all import Viewer
        from openalea.mtg.plantframe.color import colormap
        
        colors = None
        if property:
            g = colormap(g, property, cmap='jet', lognorm=True)
            colored = g.property('color')
            colors = {vid:colored.get(vid,colored.get(g.complex(vid),[0,0,0])) for vid in g.vertices(scale=g.max_scale())}
        else:
            colors = None
        s = plot3d(g, colors=colors)  # use the one of openalea.plantframe.color instead ?
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
        
    def axis_statistics(self, g):
        df_lai = self.get_exposed_areas(g, convert=True)
        axstat = None
        if not df_lai.empty:
            axstat, _ = axis_statistics(df_lai, self.domain_area, self.convUnit)
        return axstat
        
    def plot_statistics(self, axstat=None):
        pstat=None
        if axstat is not None:
            pstat = plot_statistics(axstat, self.nplants, self.domain_area)
        return pstat
        
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
        
    def get_midribs(self,g, resample = False):
        
        vids = [vid for vid in g.vertices(scale=g.max_scale() - 1) if g.label(vid).startswith('blade')]
        visible_length = g.property('visible_length')
        midribs = {vid:self.leaves.midrib(g.node(vid), resample=resample) for vid in vids if visible_length[vid] > 0}
        #
        anchor = g.property('anchor_point')
        midribs_anchor = {vid:[anchor[cid] for cid in g.components(vid) if cid in anchor] for vid in midribs}
        hins = {k:v[0][2] + midribs[k][2] for k,v in midribs_anchor.iteritems() if len(v) > 0}

        metamer = {vid:g.complex(vid) for vid in midribs}
        axe = {vid:g.complex(metamer[vid]) for vid in midribs}
        plant = {vid:g.complex(axe[vid]) for vid in midribs}
        ntop=g.property('ntop')
        
        res = [pandas.DataFrame({'vid':vid,
                                 'ntop':ntop[vid],
                                 'metamer':int(g.label(metamer[vid]).split('metamer')[1]),
                                 'axe':g.label(axe[vid]),
                                 'plant':int(g.label(plant[vid]).split('plant')[1]),
                                  'x':midribs[vid][0],
                                  'y':midribs[vid][1],
                                  'hins':hins[vid]}) for vid in hins]#hins keys are for midribs keys wich also have a geometry (anchor point
        
        return pandas.concat(res)


        
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
        