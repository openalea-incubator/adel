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

    def __init__(self, nplants=1, nsect=1, seed=None, leaves=None, stand=None, split=False,
                 face_up=False, classic=False, scene_unit='cm',
                 leaf_db=None, positions=None):

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
        self.canopy_age = 0
        self.seed = seed

        self.new_stand(nplants, seed)

    def new_stand(self, nplants=None, seed=None):
        if seed is not None:
            self.seed = seed
            numpy.random.seed(self.seed)
        if nplants is not None:
            self.nplants, self.domain, self.positions, \
            self.domain_area = self.stand.smart_stand(
                nplants, convunit=1. / self.convUnit)


    def get_axis(self, g, plant='plant1', axe='MS'):
        """ return a new mtg representing an axe
        """
        p = [vid for vid in g.vertices(scale=1) if g.label(vid) == plant][0]
        ax = [vid for vid in g.components(p) if g.label(vid) == axe][0]
        return g.sub_mtg(ax, copy=True)

    def plot(self, g, property=None):
        s = self.scene(g,property)
        Viewer.display(s)
        return s

    def scene(self, g, property=None):
        if property:
            g = colormap(g, property, cmap='jet', lognorm=True)
            colored = g.property('color')
            colors = {
            vid: colored.get(vid, colored.get(g.complex(vid), [0, 0, 0])) for
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
        pstat = None
        if axstat is not None:
            pstat = plot_statistics(axstat, self.nplants, self.domain_area)
        return pstat

    def save(self, g, index=0, dir='./adel_saved'):
        if not os.path.exists(dir):
            os.mkdir(dir)
        s = self.scene(g)
        geom = {sh.id: sh.geometry for sh in s}
        g.remove_property('geometry')
        fgeom = dir + '/scene%04d.bgeom' % (index)
        fg = dir + '/adel%04d.pckl' % (index)
        s.save(fgeom, 'BGEOM')
        f = open(fg, 'w')
        pickle.dump([g, self.canopy_age], f)
        f.close()
        # restore geometry
        g.add_property('geometry')
        g.property('geometry').update(geom)
        return fgeom, fg

    def load(self, index=0, dir='./adel_saved'):
        fgeom = dir + '/scene%04d.bgeom' % (index)
        fg = dir + '/adel%04d.pckl' % (index)

        s = Scene()
        s.read(fgeom, 'BGEOM')
        geom = {sh.id: sh.geometry for sh in s}

        f = open(fg)
        g, TT = pickle.load(f)
        f.close()

        self.canopy_age = TT

        g.add_property('geometry')
        g.property('geometry').update(geom)

        return g, TT

    def get_midribs(self, g, resample=False):

        vids = [vid for vid in g.vertices(scale=g.max_scale() - 1) if
                g.label(vid).startswith('blade')]
        visible_length = g.property('visible_length')
        midribs = {vid: self.leaves.midrib(g.node(vid), resample=resample) for
                   vid in vids if visible_length[vid] > 0}
        #
        anchor = g.property('anchor_point')
        midribs_anchor = {
        vid: [anchor[cid] for cid in g.components(vid) if cid in anchor] for vid
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
