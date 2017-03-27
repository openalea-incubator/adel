"""Mother class for Adel models"""
import os
import pandas
import numpy
import warnings

try:
    import cPickle as pickle
except:
    import pickle
from openalea.plantgl.all import Viewer, Scene
from openalea.mtg.plantframe.color import colormap

from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.mtg_interpreter import plot3d
from alinea.adel.postprocessing import axis_statistics, plot_statistics, \
    midrib_statistics
from alinea.adel.newmtg import exposed_areas, exposed_areas2canS


class Adel(object):
    """Mother class for adel models"""
    conv_units = {'mm': 0.001, 'cm': 0.01, 'dm': 0.1, 'm': 1, 'dam': 10,
                  'hm': 100,
                  'km': 1000}

    def __init__(self, nplants=1, nsect=1, seed=None, leaves=None, stand=None,
                 split=False,
                 face_up=False, classic=False, scene_unit='cm',
                 leaf_db=None, positions=None, age=0):

        if leaf_db is not None:
            warnings.warn(
                '!!!!Warning!!!! leaf_db argument is deprecated, '
                'use adel.geometric_elements.Leaves class instead',
                DeprecationWarning)
        if positions is not None:
            warnings.warn(
                '!!!!Warning!!!! positions argument is deprecated,'
                ' use stand = adel.Stand class instead',
                DeprecationWarning)

        if leaves is None:
            leaves = Leaves()

        if stand is None:
            stand = AgronomicStand(sowing_density=250, plant_density=250,
                                   inter_row=0.15)

        self.stand = stand
        self.leaves = leaves
        self.scene_unit = scene_unit
        self.convUnit = self.conv_units[self.scene_unit]
        self.nsect = nsect
        self.split = split
        self.face_up = face_up
        self.classic = classic
        self.seed = seed
        self.meta = {}

        self.new_stand(nplants, seed)
        self.new_age(age)

    def new_stand(self, nplants=None, seed=None):
        if seed is not None:
            self.seed = seed
            numpy.random.seed(self.seed)
        if nplants is not None:
            self.nplants, self.domain, self.positions, \
            self.domain_area = self.stand.smart_stand(
                nplants, convunit=1. / self.convUnit)
            self.plant_azimuths = numpy.random.random(self.nplants) * 360
            stand_parameters = {'sowing_density': self.stand.sowing_density,
                                'plant_density': self.stand.plant_density,
                                'inter_row': self.stand.inter_row,
                                'noise': self.stand.noise,
                                'density_curve_data': self.stand.density_curve_data}
            self.meta.update(
                {'stand': stand_parameters, 'nplants': self.nplants,
                 'domain': self.domain, 'domain_area': self.domain_area,
                 'nsect': self.nsect, 'scene_unit': self.scene_unit,
                 'convUnit': self.convUnit,
                 'split': self.split})

    def new_age(self, age):
        self.canopy_age = age
        self.meta.update({'canopy_age': age})

    def meta_informations(self, g):
        if 'meta' in g.property_names():
            return g.property('meta').values()[0]
        else:
            return self.meta

    @staticmethod
    def get_axis(g, plant='plant1', axe='MS'):
        """ return a new mtg representing an axe
        """
        p = [vid for vid in g.vertices(scale=1) if g.label(vid) == plant][0]
        ax = [vid for vid in g.components(p) if g.label(vid) == axe][0]
        return g.sub_mtg(ax, copy=True)

    @staticmethod
    def plot(g, property=None):
        s = Adel.scene(g, property)
        Viewer.display(s)
        return s

    @staticmethod
    def scene(g, property=None):
        if property:
            g = colormap(g, property, cmap='jet', lognorm=True)
            colored = g.property('color')
            colors = {
                vid: colored.get(vid, colored.get(g.complex(vid), [0, 0, 0]))
                for
                vid in g.vertices(scale=g.max_scale())}
        else:
            colors = None
        s = plot3d(g,
                   colors=colors)  # use the one of openalea.plantframe.color instead ?
        return s

    def get_exposed_areas(self, g, convert=False, TT=None):
        areas = exposed_areas(g)
        if convert:
            areas = exposed_areas2canS(areas)
        if TT is None:
            TT = self.meta_informations(g)['canopy_age']
        areas['TT'] = TT
        return areas

    def axis_statistics(self, g):
        meta = self.meta_informations(g)
        df_lai = self.get_exposed_areas(g, convert=True)
        axstat = None
        if not df_lai.empty:
            axstat, _ = axis_statistics(df_lai, meta['domain_area'],
                                        meta['convUnit'])
        return axstat

    def plot_statistics(self, g=None, axstat=None):
        meta = self.meta_informations(g)
        pstat = None
        if axstat is None:
            axstat = self.axis_statistics(g)
        pstat = plot_statistics(axstat, meta['nplants'], meta['domain_area'])
        return pstat

    def save(self, g, index=0, dir='./adel_saved', basename=None,
             check_meta=True):
        if check_meta:
            if 'meta' not in g.property_names():
                root = g.node(0)
                root.meta = self.meta
        if basename is None:
            if not os.path.exists(dir):
                os.mkdir(dir)
            basename_geom = dir + '/scene%04d' % (index)
            basename_adel = dir + '/adel%04d' % (index)
        else:
            basename_adel = basename_geom = str(basename)
        s = Adel.scene(g)
        geom = {sh.id: sh.geometry for sh in s}
        g.remove_property('geometry')
        fgeom = basename_geom + '.bgeom'
        fg = basename_adel + '.pckl'
        s.save(fgeom, 'BGEOM')
        with open(fg, 'w') as output:
            pickle.dump(g, output)
        # restore geometry
        g.add_property('geometry')
        g.property('geometry').update(geom)
        return fgeom, fg

    @staticmethod
    def load(index=0, dir='./adel_saved', basename=None, load_geom=True):
        if basename is None:
            if not os.path.exists(dir):
                os.mkdir(dir)
            basename_geom = dir + '/scene%04d' % (index)
            basename_adel = dir + '/adel%04d' % (index)
        else:
            basename_adel = basename_geom = basename
        fgeom = basename_geom + '.bgeom'
        fg = basename_adel + '.pckl'
        if not os.path.exists(fgeom) or not os.path.exists(fg):
            raise IOError('adel cannot find saved files')

        f = open(fg)
        g = pickle.load(f)
        f.close()

        if load_geom:
            s = Scene()
            s.read(fgeom, 'BGEOM')
            geom = {sh.id: sh.geometry for sh in s}
            g.add_property('geometry')
            g.property('geometry').update(geom)

        return g

    def get_midribs(self, g, resample=False):

        vids = [vid for vid in g.vertices(scale=g.max_scale() - 1) if
                g.label(vid).startswith('blade')]
        visible_length = g.property('visible_length')
        midribs = {vid: self.leaves.midrib(g.node(vid), resample=resample) for
                   vid in vids if visible_length[vid] > 0}
        #
        anchor = g.property('anchor_point')
        midribs_anchor = {
            vid: [anchor[cid] for cid in g.components(vid) if cid in anchor] for
            vid
            in midribs}
        hins = {k: v[0][2] + midribs[k][2] for k, v in
                midribs_anchor.iteritems() if len(v) > 0}

        metamer = {vid: g.complex(vid) for vid in midribs}
        axe = {vid: g.complex(metamer[vid]) for vid in midribs}
        plant = {vid: g.complex(axe[vid]) for vid in midribs}
        ntop = g.property('ntop')

        res = [pandas.DataFrame({'vid': vid,
                                 'ntop': ntop[vid],
                                 'metamer': int(
                                     g.label(metamer[vid]).split('metamer')[1]),
                                 'axe': g.label(axe[vid]),
                                 'plant': int(
                                     g.label(plant[vid]).split('plant')[1]),
                                 'x': midribs[vid][0],
                                 'y': midribs[vid][1],
                                 'hins': hins[vid]}) for vid in
               hins]  # hins keys are for midribs keys wich also have a geometry (anchor point

        return pandas.concat(res)

    def midrib_statistics(self, g):
        data = self.get_midribs(g)
        return midrib_statistics(data)
