""" Class instanciating a wheat canopy and complying to astk canopy interface
"""

import pandas

from alinea.adel.AdelR import setAdel, RunAdel, genGeoAxe, checkAxeDyn, getAxeT, \
    getPhenT, getPhytoT, saveRData
from alinea.adel.newmtg import move_properties
import alinea.adel.data_samples as adel_data
from alinea.adel.mtg_interpreter import plot3d
from alinea.astk.TimeControl import *


from alinea.adel.adel import Adel


class AdelWheat(Adel):

    def __init__(self, devT=None, sample='random', thermal_time_model=None,
                 geoAxe=None, incT=60, dinT=5, dep=7, run_adel_pars=None,
                 aborting_tiller_reduction=1.0, ssipars=None,
                 nplants=1, duplicate=None, species=None, nsect=1, leaves=None, stand=None,
                 aspect='smart', split=False,
                 face_up=False, classic=False, devT_unit='cm', scene_unit='cm',
                 age=None, seed=None,
                 leaf_db=None,
                 positions=None,
                 convUnit=None):

        if species is not None or isinstance(leaves, dict):
            raise ValueError('multi_species canopies not yet implemented')

        if devT is None:
            devT = adel_data.devT()
            devT_unit = 'cm'

        if devT_unit != scene_unit:
            convert = self.conv_units[devT_unit] / self.conv_units[scene_unit]
            for x in (
            'L_internode', 'W_internode', 'W_sheath', 'W_blade', 'L_sheath',
            'L_blade'):
                devT['dimT'][x] *= convert

        self.devT = devT

        self.ref_plants = list(set(self.devT['axeT']['id_plt']))
        super(AdelWheat, self).__init__(nref_plants=len(self.ref_plants),
                                        nplants=nplants, duplicate=duplicate,
                                        nsect=nsect, species=species,
                                        leaves=leaves, stand=stand,
                                        aspect=aspect,
                                        split=split, face_up=face_up,
                                        classic=classic, scene_unit=scene_unit,
                                        age=age,
                                        seed=seed, leaf_db=leaf_db,
                                        positions=positions,
                                        convUnit=convUnit)


        if run_adel_pars is None:
            run_adel_pars = {'senescence_leaf_shrink': 0.5, 'leafDuration': 2,
                             'fracLeaf': 0.2, 'stemDuration': 2. / 1.2,
                             'dHS_col': 0.2, 'dHS_en': 0, 'epsillon': 1e-6,
                             'HSstart_inclination_tiller': 1,
                             'rate_inclination_tiller': 30, 'drop_empty': True}

        if thermal_time_model is None:
            thermal_time_model = DegreeDayModel(Tbase=0)

        if geoAxe is None:
            geoAxe = genGeoAxe(incT=incT, dinT=dinT, dep=dep)

        assert len(self.leaves.keys()) == 1
        k = self.leaves.keys()[0]

        pars = {'devT': self.devT, 'RcodegeoLeaf': self.leaves[k].geoLeaf,
                    'RcodegeoAxe': geoAxe,
                    'seed': self.seed, 'sample': sample, 'xydb': self.leaves[k].xydb,
                    'srdb': self.leaves[k].srdb, 'ssipars': ssipars}

        if self.duplicate is None:
            self.pars = setAdel(nplants=self.nplants, **pars)
        else:
            if self.nquot > 0:
                self.pars_quot = setAdel(nplants=self.nquot, **pars)
            if self.nrem > 0:
                self.pars_rem = setAdel(nplants=self.nrem, **pars)

        self.thermal_time = thermal_time_model
        self.run_adel_pars = run_adel_pars
        self.aborting_tiller_reduction = aborting_tiller_reduction

    @staticmethod
    def timing(delay, steps, weather, start_date):
        """ compute timing and time_control_sets for a simulation between start and stop. 

        Return 0 when there is no rain
        """
        timestep = delay
        start_date = weather.str_to_datetime(start_date)
        temp = [d['temperature_air'] for d in
                weather.split_weather(timestep, start_date, steps)]

        return (TimeControlSet(Tair=temp[i / int(delay)],
                               dt=delay) if not i % delay  else TimeControlSet(
            dt=0) for i in range(steps))

    def setup_canopy(self, age=None):
        """ Create a canopy of a given age

        Args:
            age: thermal time age of the canopy (simce emergence 1).
            If None (default) a full grown canopy is returned (TT = last
            ligulation + 300)

        Returns:

        """

        if age is None:
            axT = pandas.DataFrame(self.devT['axeT'])
            phenT = pandas.DataFrame(self.devT['phenT'])
            axT = axT.loc[axT['TT_stop_axis'] == 'NA']
            axT = axT.merge(phenT)
            TT_last_col = (axT['TT_col_phytomer1'] + axT['dTT_col_phytomer']).max()
            age = TT_last_col + 300


        self.new_stand(age=age)

        if self.duplicate is None:
            canopy = RunAdel(age, self.pars, adelpars=self.run_adel_pars)
            # add index_relative_to_MS_phytomer to parameters
            if 'index_relative_to_MS_phytomer' in self.devT['dimT']:
                dfd = pandas.DataFrame(self.devT['axeT']).merge(
                    pandas.DataFrame(self.devT['dimT'])).loc[:, (
                      'id_plt', 'id_axis', 'index_phytomer',
                      'index_relative_to_MS_phytomer')]
                canopy = pandas.DataFrame(canopy).merge(dfd.rename(
                    columns={'id_plt': 'plant', 'id_axis': 'axe_id',
                             'index_phytomer': 'numphy'})).to_dict('list')
            stand = zip(self.positions, self.plant_azimuths)
            g = self.build_mtg(canopy, stand,
                               aborting_tiller_reduction=self.aborting_tiller_reduction)
        else:
            # produce plants positionned at origin
            grem = None
            if self.nrem > 0:
                canopy = RunAdel(age, self.pars_rem, adelpars=self.run_adel_pars)
                grem = self.build_mtg(canopy, stand=None,
                               aborting_tiller_reduction=self.aborting_tiller_reduction)

            if self.nquot > 0:
                canopy = RunAdel(age, self.pars_quot, adelpars=self.run_adel_pars)
                gquot = self.build_mtg(canopy, stand=None,
                               aborting_tiller_reduction=self.aborting_tiller_reduction)

            g = self.duplicated(gquot, grem)

        return g

    def checkAxeDyn(self, dates=range(0, 2000, 100), density=None):
        if density is None:
            density = self.stand.plant_density
        return checkAxeDyn(self.pars, dates, density)

    def check_nff(self):
        """ Count the probability of occurence of MS with given nff
        """
        ms = numpy.array(self.devT['axeT']['id_axis']) == 'MS'
        nffs = numpy.array(self.devT['axeT']['N_phytomer'])[ms]
        counts = {n: nffs.tolist().count(n) for n in set(nffs)}
        probas = {n: nffs.tolist().count(n) * 1.0 / len(nffs) for n in
                  set(nffs)}
        return counts, probas

    def check_primary_tillers(self):
        """ Count/estimate probabilitie of occurence of primary tillers
        """
        import re
        axis = self.devT['axeT']['id_axis']
        ms = [e for e in axis if re.match('MS', e)]
        tillers = [e for e in axis if re.match('T.$', e)]
        counts = {n: tillers.count(n) for n in set(tillers)}
        probas = {n: tillers.count(n) * 1.0 / len(ms) for n in set(tillers)}
        return counts, probas

    def axeT(self):
        return getAxeT(self.pars)

    def phenT(self, axe='MS'):
        return getPhenT(self.pars, axe=axe)

    def phytoT(self, axe='MS'):
        return getPhytoT(self.pars, axe=axe)

    def save_pars(self):
        saveRData(self.pars, 'plants', 'adel_pars.RData')

    def grow(self, g, time_control):

        try:  # old interface
            if time_control.dt <= 0:
                dday = 0.
            else:
                dday = time_control.Tair.mean()
        except:
            data = time_control
            tt = self.thermal_time(data.index, data)
            dday = tt[-1]

        # refg = self.setup_canopy(age = self.canopy_age)
        self.canopy_age += dday
        newg = self.setup_canopy(age=self.canopy_age)
        # newg = mtg_update(newg, g, refg)
        move_properties(g, newg)

        return newg

    def grow_dd(self, g, dday):
        # refg = self.setup_canopy(age = self.canopy_age)
        self.canopy_age += dday
        newg = self.setup_canopy(age=self.canopy_age)
        # newg = mtg_update(newg, g, refg)
        move_properties(g, newg)
        return newg


def adelwheat_node(nplants=1, nsect=1, devT=None, leaves=None, geoAxe=None,
                   stand=None, run_adel_pars=None, options={}):
    args = locals()
    args.update(options)
    args.pop('options')
    model = AdelWheat(**args)
    return model


def setup_canopy_node(adel, age=10):
    return adel.setup_canopy(age)


def adel_scene_node(g):
    return plot3d(g)


def axis_statistics_node(adel, g):
    return adel.axis_statistics(g),


def plot_statistics_node(adel, axstat):
    return adel.plot_statistics(axstat),


# user friendly macros
from alinea.adel.stand.stand import agronomicplot
from alinea.astk.plant_interface import *


def initialise_stand(age=0., length=0.1, width=0.2, sowing_density=150,
                     plant_density=150, inter_row=0.12, nsect=1, seed=None,
                     sample='random'):
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
    nplants, positions, domain, domain_area, convUnit = agronomicplot(
        length=length,
        width=width,
        sowing_density=sowing_density,
        plant_density=plant_density,
        inter_row=inter_row)
    wheat = AdelWheat(nplants=nplants, positions=positions, nsect=nsect,
                      seed=seed, sample=sample)
    g, _ = new_canopy(wheat, age=age)
    return g, wheat, domain_area, domain, convUnit
