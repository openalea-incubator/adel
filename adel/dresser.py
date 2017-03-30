"""User interface for static reconstruction using Adel modules
"""
import numpy
import pandas

from alinea.adel.geometric_elements import Leaves
from alinea.adel.newmtg import mtg_factory
from alinea.adel.mtg_interpreter import mtg_interpreter
from alinea.adel.adel import Adel


def blade_dimension(area=None, length=None, width=None, ntop=None, leaves=None,
                    plant=1, wl=0.1):
    """Estimate blade dimension and/or compatibility with leaf shapes form factors

    Args:
        area: (array) vector of blade area. If None, will be estimated using
         other args
        length: (array) vector of blade lengths. If None, will be estimated using
         other args
        width: (array) vector of blade widths. If None, will be estimated using
         other args
        ntop: (array) vector of leaf position (topmost leaf =1). If None
        (default), leaf dimensions are assumed to be from top to base.
        leaves: (object) a Leaves instance defining leaf shapes. If None (default)
        the adel default shape will be used
        plant: (int or array) vector of plant number
        wl: (float) the width / length ratio used to estimates dimensions in case
         of  uncomplete data

    Returns:
        a pandas dataframe with estimated blade dimensions

    """

    if area == length == width == None:
        area = (15, 20, 30)

    if leaves is None:
        leaves = Leaves()
    ff = leaves.form_factor()
    ff = {int(k): v for k, v in ff.iteritems()}

    if area is None:
        if length is None:
            width = numpy.array(width)
            length = width / numpy.array(wl)
        elif width is None:
            length = numpy.array(length)
            width = length * numpy.array(wl)
        else:
            length = numpy.array(length)
            width = numpy.array(width)
        if ntop is None:
            ntop = numpy.arange(1, len(length) + 1)
        else:
            ntop = numpy.array(ntop)
        ffn = numpy.array([ff[k] for k in ntop])
        area = ffn * length * width
    else:
        area = numpy.array(area)
        if ntop is None:
            ntop = numpy.arange(1, len(area) + 1)
        else:
            ntop = numpy.array(ntop)
        ffn = numpy.array([ff[k] for k in ntop])
        # adjust length/width if one is  None or overwrite width if all are set
        if length is None:
            if width is None:
                length = numpy.sqrt(area / ffn / wl)
                width = length * wl
            else:
                width = numpy.array(width)
                length = area / ffn / width
        else:
            length = numpy.array(length)
            width = area / ffn / length

    if isinstance(plant, int):
        plant = [plant] * len(ntop)

    return pandas.DataFrame({'plant': plant, 'ntop': ntop, 'L_blade': length,
                             'W_blade': width, 'S_blade': area})


def stem_dimension(h_ins=None, d_stem=None, internode=None, sheath=None,
                   d_internode=None, d_sheath=None, ntop=None, plant=1):
    """Estimate dimension of stem organs from stem measurements

    Args:
        h_ins: (array) vector of blade insertions height
        d_stem:(float or array) vector of stem diameter
        internode:(array) vector of internode lengths. If None, will be estimated using
         other args
        sheath: (array) vector of sheath lengths. If None, will be estimated using
         other args
        d_internode: (array) vector of intenode diameters. If None, will be estimated using
         other args
        d_sheath: (array) vector of sheath diameters. If None, will be estimated using
         other args
        ntop:(array) vector of leaf position (topmost leaf =1). If None
        (default), stem dimensions are assumed to be from top to base.
        plant: (int or array) vector of plant number

    Returns:
        a pandas dataframe with estimated sheath and internode dimension
    """

    if h_ins is None and h_ins == internode == sheath:
        h_ins = (60, 50, 40)

    if d_stem is None and d_stem == d_internode == d_sheath:
        d_stem = 0.3

    if h_ins is None:
        if sheath is None:
            sheath = numpy.array([0] * len(internode))
        else:
            sheath = numpy.array(sheath)
        if internode is None:
            internode = numpy.array([0] * len(sheath))
        else:
            internode = numpy.array(internode)
        if ntop is None:
            ntop = numpy.arange(1, len(h_ins) + 1)
        else:
            ntop = numpy.array(ntop)
        order = numpy.argsort(-ntop)
        reorder = numpy.argsort(order)
        h_ins = (internode[order].cumsum() + sheath[order])[reorder]
    else:
        h_ins = numpy.array(h_ins)
        if ntop is None:
            ntop = numpy.arange(1, len(h_ins) + 1)
        else:
            ntop = numpy.array(ntop)
        order = numpy.argsort(-ntop)
        reorder = numpy.argsort(order)

        if sheath is None:
            if internode is None:
                sheath = numpy.array([0] * len(h_ins))
                internode = numpy.diff([0] + list(h_ins[order]))[reorder]
            else:
                internode = numpy.array(internode)
                sheath = \
                numpy.maximum(0, h_ins[order] - internode[order].cumsum())[
                    reorder]
                internode = \
                numpy.diff([0] + (h_ins[order] - sheath[order]).tolist())[
                    reorder]
        else:
            sheath = numpy.array(sheath)
            internode = \
            numpy.diff([0] + (h_ins[order] - sheath[order]).tolist())[reorder]

    if d_internode is None:
        if d_sheath is None:
            d_internode = [d_stem] * len(h_ins)
        else:
            d_internode = d_sheath
    if d_sheath is None:
        d_sheath = d_internode

    if isinstance(plant, int):
        plant = [plant] * len(ntop)

    return pandas.DataFrame({'plant': plant, 'ntop': ntop, 'h_ins': h_ins,
                             'L_sheath': sheath, 'W_sheath': d_sheath,
                             'L_internode': internode,
                             'W_internode': d_internode})


def ear_dimension(peduncle=None, ear=None, spike=None, d_peduncle=0.3,
                  projected_area_ear=None, wl_ear=0.1, plant=1):
    """Estimate dimensions of ear cylinders from ear measuremenst

    Args:
        peduncle: length of the peduncle. If None, no peduncle is computed
        ear: length of the ear. If None, no ear is computed
        spike: length of the ear + awn. If None, no awn is computed
        d_peduncle: diameter of the peduncle
        projected_area_ear: projected area of the ear
        wl_ear: width/length ratio for ear. Used only if projected_ear_area is missing
        plant: (int or array) vector of plant number

    Returns:
        a pandas DataFrame with ear dimensions
    """
    if peduncle == ear == spike == None:
        return None

    dfl = []
    pos = 0

    if peduncle is not None:
        dfl.append(pandas.DataFrame({'plant': plant, 'ntop': pos,
                                     'L_internode': peduncle,
                                     'W_internode': d_peduncle}, index=[pos]))
        pos -= 1

    if ear is not None:
        if projected_area_ear is None:
            w_ear = wl_ear * ear
        else:
            w_ear = float(projected_area_ear) / ear
        dfl.append(pandas.DataFrame(
            {'plant': plant, 'ntop': pos, 'L_internode': ear,
             'W_internode': w_ear}, index=[pos]))
        pos -= 1

        if spike is not None:
            dfl.append(pandas.DataFrame(
                {'plant': plant, 'ntop': pos, 'L_internode': spike - ear,
                 'W_internode': w_ear}, index=[pos]))

    return pandas.concat(dfl, axis=0)


def dimension_table(blades=None, stem=None, ear=None):
    if blades is None:
        blades = blade_dimension()
    if stem is None:
        stem = stem_dimension()

    if ear is None:
        return blades.merge(stem)
    else:
        stemear = pandas.concat([stem, ear]).set_index(['plant', 'ntop'])
        return pandas.concat([stemear, blades.set_index(['plant', 'ntop'])],
                             axis=1).reset_index()


class AdelDress(Adel):
    """A class interface to Adel for static reconstruction"""

    def __init__(self, dimT=None, dim_unit='cm', nplants=1, duplicate=None,
                 nsect=1,
                 leaves=None, stand=None,
                 aspect='smart', split=False,
                 face_up=False, classic=False, scene_unit='cm', age=None,
                 seed=None):
        """ Instantiate a dresser

        Args:
            dimT: (panda.dataFrame) a table with organ dimensions
            dim_unit: (string) length unit used in the dimension table
            nsect: (int) the number of sectors on leaves
            leaves: (object) a Leaves class instance pointing to leaf shape
             database or a {species:leaf_db} dict referencing distinct database
             per species
            stand: (object) a Stand class instance
            aspect: (str) the aspect of the stand (square, line or smart)
            classic: (bool) should stem cylinders be classical pgl cylinders ?
            scene_unit: (string) desired length unit for the output mtg
            age: (optional) the age of the canopy
        """
        super(AdelDress, self).__init__(nplants=nplants, duplicate=duplicate,
                                        nsect=nsect,
                                        leaves=leaves, stand=stand,
                                        aspect=aspect,
                                        split=split, face_up=face_up,
                                        classic=classic, scene_unit=scene_unit,
                                        age=age,
                                        seed=seed)

        if dimT is None:
            dimT = dimension_table()

        self.dim_unit = dim_unit
        self.dimT = dimT.fillna(0)
        self.dimT = dimT

    def canopy_table(self, nplants=None, azimuth=None, species=None):

        if nplants is None:
            nplants = self.nplants

        if azimuth is None:
            def azimuth(n, ntop, axe):
                return 180 + (numpy.random.random() - 0.5) * 30

        if species is None:
            species = {0: 1}

        df = self.dimT.loc[:,
             ['plant', 'ntop', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath',
              'L_internode', 'W_internode']]
        df.rename(
            columns={'L_blade': 'Ll', 'W_blade': 'Lw_shape',
                     'L_sheath': 'Gl', 'W_sheath': 'Gd',
                     'L_internode': 'El', 'W_internode': 'Ed'}, inplace=True)
        conv = self.conv_units[self.dim_unit] / self.conv_units[self.scene_unit]
        df.loc[:, ('Ll', 'Lw_shape', 'Gl', 'Gd', 'El', 'Ed')] *= conv
        # add mandatory topological info and sort from base to top
        df.loc[:, 'axe_id'] = 'MS'
        df.loc[:, 'ms_insertion'] = 0
        df.loc[:, 'numphy'] = df.ntop.max() + 1 - df.ntop
        # additions for statistics
        df.loc[:, 'nff'] = df['numphy'].max()
        df.loc[:, 'HS_final'] = df['numphy'].max()
        df = df.sort_values(['plant', 'numphy'])

        ref_plant = numpy.random.choice(list(set(df['plant'])), nplants)
        labels = species.keys()
        nbspec = len(labels)
        spec = [labels[k] for k in
         numpy.random.choice(nbspec, nplants, p=species.values())]
        dfl = []
        for i, p in enumerate(ref_plant):
            dfp = df.loc[df['plant'] == p, :]
            dfp['refplant_id'] = p
            dfp['species'] = spec[i]
            dfp.loc[:, 'plant'] = i + 1
            # compute visibility
            ht0 = 0
            hbase = dfp['El'].cumsum() - dfp['El']
            hcol = hbase + dfp['Gl'] + dfp['El']
            h_hide = [max([ht0] + hcol[:i].tolist()) for i in range(len(hcol))]
            htube = numpy.maximum(0, h_hide - hbase)
            dfp['Lv'] = numpy.minimum(dfp['Ll'], numpy.maximum(0,
                                                               dfp['Ll'] + dfp[
                                                                   'Gl'] + dfp[
                                                                   'El'] - htube))
            dfp['Gv'] = numpy.minimum(dfp['Gl'], numpy.maximum(0,
                                                               dfp['Gl'] + dfp[
                                                                   'El'] - htube))
            dfp['Ev'] = numpy.minimum(dfp['El'],
                                      numpy.maximum(0, dfp['El'] - htube))
            # add missing mandatory data  (does like adel)
            dfp.loc[:, 'Laz'] = [azimuth(*arg) for arg in
                                 zip(dfp['numphy'], dfp['ntop'],
                                     dfp['axe_id'])]  # leaf azimuth
            dfp.loc[:, 'LcType'] = numpy.where(dfp['ntop'] > 0, dfp['ntop'],
                                               1)  # selector for first level in leaf db
            dfp.loc[:,
            'LcIndex'] = 1 + numpy.array(map(lambda (s,t): numpy.random.choice(
                range(len(self.leaves[s].xydb[str(t)]))), zip(dfp['species'], dfp[
                                                 'LcType'])))  # selector for second level in leaf_db (ranging 1:max_nb_leaf_per_level)
            # fill other columns
            dfp.loc[:, 'Lr'] = 0
            dfp.loc[:, 'Lsen'] = 0
            dfp.loc[:, 'L_shape'] = dfp['Ll']
            dfp.loc[:, 'Linc'] = 1
            dfp.loc[:, 'Gsen'] = 0
            dfp.loc[:, 'Ginc'] = 0
            dfp.loc[:, 'Esen'] = 0
            dfp.loc[:, 'Einc'] = 0
            dfl.append(dfp)

        return pandas.concat(dfl)

    def canopy(self, nplants=None, duplicate=None, azimuth=None, species=None,
               seed=None,
               age=None, aspect='smart'):
        """ Generate a mtg encoding the canopy

        Args:
            nplants: the number of plants in the canopy
            duplicate:
            azimuth: a callable returning leaf azimuth as a function of leaf rank,
             leaf rank from top and axe_id. If None
            species : a {species: frequency} dict indicating the composition of
            the canopy. If None, a monospecific canopy of species '0' is generated
            seed: (int) a value to initialize random number generator

        Returns:

        """

        self.new_stand(nplants=nplants, duplicate=duplicate, seed=seed,
                       aspect=aspect, age=age)
        df = self.canopy_table(azimuth=azimuth, species=species)

        stand = zip(self.positions, self.plant_azimuths)

        if self.duplicate is None:
            g = mtg_factory(df.to_dict('list'), leaf_sectors=self.nsect,
                            leaves=self.leaves, stand=stand)
            g = mtg_interpreter(g, self.leaves, classic=self.classic,
                                face_up=self.face_up)
        else:
            raise NotImplementedError("duplication not yet implemented for dresser")

        return g
