# -*- python -*-
#
#       Copyright 2015 INRIA - CIRAD - INRA
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       WebSite : https://github.com/openalea-incubator/adel
#
# ==============================================================================
"""
A place to develop candidates methods for mtg.py / topological builder

"""
from alinea.adel.exception import *

# temporary import
from alinea.adel.mtg import convert, properties_from_dict

# import csv

from openalea.mtg import MTG, fat_mtg
from openalea.mtg.traversal import iter_mtg

import numpy
import pandas


# to do add hypocotyl
# try:
# from openalea.mtg.traversal import *
# except:
# pass

# from openalea.mtg.io import read_lsystem_string
# from openalea.mtg.algo import union

# import numpy
# try:
# from openalea.plantgl.all import (Scene,Translation, Vector3, Geometry,
# AxisRotation, AxisRotated,  Transform4, BaseOrientation,
# Shape, Material, Color3, PglTurtle, Mesh, Translated)
# except:
# Material = tuple
# Color3 = tuple
# pass
# import random
# from math import pi

def internode_elements(l, lvis, lsen, az, inc, d, split=False):
    """ returns parameters of internode elements (visible parts of the internode).
    l is the length of the internode
    lv is the visible length (senesced + green)
    lsen is the senescent apical length
    az is the azimuth angle
    inc is the inclination angle
    Fisrt element is for the green visible part of the internode
    Second element is for senesced visible part of the internode
    """
    lhide = None
    lgreen = None
    is_green = None
    Svis = None
    Sgreen = None
    Ssen = None
    Shide = None
    try:
        lhide = max(l - lvis, 0.)
        lgreen = lvis - min(lsen, lvis)
        lsen = lvis - lgreen
        #
        Shide = numpy.pi * lhide * d
        Svis = numpy.pi * lvis * d
        Sgreen = numpy.pi * lgreen * d
        Ssen = numpy.pi * lsen * d
        #
        is_green = lgreen >= lsen
    except TypeError:
        pass

    hidden_elt = {'label': 'HiddenElement', 'length': lhide, 'area': Shide,
                  'is_green': True}
    if split:
        green_elt = {'label': 'StemElementg', 'length': lgreen,
                     'is_green': True, 'azimuth': az, 'inclination': inc}
        sen_elt = {'label': 'StemElements', 'length': lsen, 'is_green': False,
                   'azimuth': 0, 'inclination': 0}
        return [hidden_elt, green_elt, sen_elt]
    else:
        elt = {'label': 'StemElement', 'length': lvis, 'area': Svis,
               'green_length': lgreen, 'green_area': Sgreen,
               'senesced_length': lsen, 'senesced_area': Ssen,
               'is_green': is_green, 'azimuth': az, 'inclination': inc}
        return [hidden_elt, elt]


def sheath_elements(l, lvis, lsen, az, inc, d, split=False):
    """ returns parameters of sheath elements (visible parts of the sheath).
    l is the length of the sheath
    lv is the visible length (senesced + green)
    lsen is the senescent apical length
    Fisrt elements is for green visible part of the sheath
    Second elements is for senesced visible part of the sheath
    """
    # same logic as internodes
    return internode_elements(l, lvis, lsen, az, inc, d, split)


def blade_elements(sectors, l, lvis, lrolled, lsen, Lshape, Lwshape, shape_key,
                   leaves, split=False):
    """ return parameters of blade elements (visible parts of the blade).
    sectors is the number of sectors dividing pattern blade shape
    l is the length of the blade
    lvis is the visible length (senesced + green, rolled + flat)
    lrolled is the visible rolled length of the blade
    lsen is the senescent apical length
    Lshape is length of the blade used as a pattern shape
    """
    lrolled = max(0, min(lrolled, lvis))
    if lrolled < 1e-6:
        lrolled = 0
    lhide = None
    lgreen = None
    lflat = None
    # s (on mature shape) at which leaf becomes flat and visible
    s_limvis = Lshape
    # s (on mature shape) at which leaf becomes senescent
    s_limsen = Lshape
    # s(on mature shape) at which leaf becomes rolled
    s_limrolled = Lshape
    # compute partitioning of visible length
    try:
        lhide = max(l - lvis, 0.)
        S_hide = 0
        if lhide > 0:
            S_hide = leaves.blade_elt_area(shape_key, Lshape, Lwshape,
                                           (Lshape - l) / Lshape,
                                           (Lshape - lvis) / Lshape)
        lflat = lvis - min(lrolled, lvis)
        lrolled = lvis - lflat
        lgreen = lvis - min(lsen, lvis)
        lsen = lvis - lgreen
        s_limvis = Lshape - lvis
        s_limsen = Lshape - lsen
        s_limrolled = Lshape - lflat
        # print(lgreen,lsen,s_limvis,s_limsen)
    except TypeError:
        pass

    # hidden part + rolled : area are set to zero as leaf rolled area is already counted in leaf element 
    # hidden_elt = {'label': 'StemElement', 'offset': lhide, 'length': lrolled, 'area': 0, 'green_length': lrolled, 'green_area': 0, 'senesced_length' : 0, 'senesced_area':0,'is_green': True, 'azimuth': 0, 'inclination':0}
    hidden_elt = {'label': 'HiddenElement', 'length': lhide, 'area': S_hide,
                  'is_green': True}
    elements = []
    ds = 0
    if Lshape is not None:
        ds = float(Lshape) / sectors
    # compute sectors from leafshape bes to leaf shape top
    if lvis > 1e-6:
        elements = [hidden_elt]
        for isect in range(sectors):
            ls_vis = 0
            ls_rolled = 0
            ls_rolled_green = 0
            ls_rolled_sen = 0
            d_rolled = 0
            ls_green = 0
            ls_sen = 0
            S_green = 0
            S_sen = 0
            S_tot = 0
            srb_green = None
            srt_green = None
            srb_sen = None
            srt_sen = None

            try:
                st = (isect + 1) * ds
                ls_vis = min(ds, max(0., st - s_limvis))
                if ls_vis > 0:
                    sb = st - ls_vis
                    sb_green = sb
                    st_green = min(st, max(sb_green, s_limsen))
                    sb_sen = st_green
                    st_sen = st
                    ls_green = st_green - sb_green
                    ls_sen = st_sen - sb_sen
                    # print(sb_green,st_green,st_sen)
                    #
                    # Compute area of elements
                    if ls_green > 0:
                        S_green = leaves.blade_elt_area(shape_key, Lshape,
                                                        Lwshape,
                                                        sb_green / Lshape,
                                                        st_green / Lshape)
                    if ls_sen > 0:
                        S_sen = leaves.blade_elt_area(shape_key, Lshape,
                                                      Lwshape, sb_sen / Lshape,
                                                      st_sen / Lshape)
                    # made intergration again for avoiding fluctuations
                    # S_tot = blade_elt_area(xysr_shape, Lshape, Lwshape, sb / Lshape, st / Lshape)
                    S_tot = S_green + S_sen
                    #
                    # compute position of flat parts of the element
                    ls_flat = min(ls_vis, max(0., st - s_limrolled))
                    if ls_flat > 0:
                        sb = st - ls_flat
                        sb_green = sb
                        st_green = min(st, max(sb_green, s_limsen))
                        sb_sen = st_green
                        st_sen = st
                        srb_green = (float(sb_green) - s_limvis) / lvis
                        srt_green = (float(st_green) - s_limvis) / lvis
                        srb_sen = (float(sb_sen) - s_limvis) / lvis
                        srt_sen = (float(st_sen) - s_limvis) / lvis
                    # print(srb_green,srt_green,srb_sen,srt_sen)
                    if ls_flat < ls_vis:
                        ls_rolled = ls_vis - ls_flat
                        ls_rolled_green = min(ls_rolled, ls_green)
                        ls_rolled_sen = ls_rolled - ls_rolled_green
                        S_rolled = S_tot * ls_rolled / ls_vis
                        d_rolled = float(S_rolled) / numpy.pi / ls_rolled

            except TypeError:  # input is None
                #    print "passing"
                pass
            green_elt = {'label': 'LeafElement' + str(isect + 1) + 'g',
                         'length': ls_green, 'area': S_green, 'is_green': True,
                         'srb': srb_green, 'srt': srt_green,
                         'lrolled': ls_rolled_green, 'd_rolled': d_rolled}
            sen_elt = {'label': 'LeafElement' + str(isect + 1) + 's',
                       'length': ls_sen, 'area': S_sen, 'is_green': False,
                       'srb': srb_sen, 'srt': srt_sen, 'lrolled': ls_rolled_sen,
                       'd_rolled': d_rolled}
            elt = {'label': 'LeafElement' + str(isect + 1),
                   'length': ls_sen + ls_green, 'area': S_tot,
                   'green_length': ls_green, 'green_area': S_green,
                   'senesced_length': ls_sen, 'senesced_area': S_sen,
                   'is_green': (ls_green > ls_sen), 'srb': srb_green,
                   'srt': srt_sen, 'lrolled': ls_rolled, 'd_rolled': d_rolled}
            if split:
                elements.extend([green_elt, sen_elt])
            else:
                elements.extend([elt])
    return elements


def adel_metamer(Ll=None, Lv=None, Lr=None, Lsen=None, L_shape=None,
                 Lw_shape=None, shape_key=None, Linc=None, Laz=None, Lsect=1,
                 Gl=None, Gv=None, Gsen=None, Gd=None, Ginc=None, El=None,
                 Ev=None, Esen=None, Ed=None, Einc=None, elongation=None,
                 ntop=None, leaves=None, **kwargs):
    """ Contructs metamer elements for adel from parameters describing a static state.
    Parameters are : 
       - Ll : length of the blade
       - Lv : visible (emerged) length of blade (green + senesced, rolled + unrolled)
       - Lr : rolled part of the blade
       - Lsen : length of the senescent part of the blade (hidden + visible)       
       - L_shape : Mature length of the blade used to compute blade shape
       - Lw_shape : Maximal width of the blade used to compute blade shape
       - xysr_shape : a (x,y,s,r) tuple describing blade geometry
       - Linc : relative inclination of the base of leaf blade at the top of leaf sheath (deg) 
       - Laz : relative azimuth of the leaf
       - Lsect : the number of sectors per leaf blade
       - Gl : length of the sheath (hidden + visible)
       - Gv : emerged length of the sheath
       - Gsen : senescent length of the sheath (hidden + visible)
       - Gd : apparent diameter of the sheaths
       - Ginc : relative inclination of the sheath
       - El: length of the internode (hidden + visible)
       - Ev: emerged length of the internode 
       - Esen: senescent length of the internode (hidden + visible)
       - Ed: diameter of the internode
       - Einc : relative inclination of the internode
 
    """
    # to do add diameter and Lrolled to blade
    Eaz = Laz
    Gaz = 0
    split = kwargs.get('split', False)
    exposition = kwargs.get('exposition', 'NA')
    lifetime = kwargs.get('lifetime', 'NA')
    mtype = kwargs.get('m_type', 'vegetative')
    age = kwargs.get('age', 'NA')
    is_ligulated = kwargs.get('is_ligulated', 'NA')

    if mtype != 'vegetative':
        modules = [
            {'label': mtype, 'ntop': ntop, 'length': El, 'visible_length': Ev,
             'senesced_length': Esen, 'diameter': Ed, 'azimuth': Eaz,
             'inclination': Einc,
             'elements': internode_elements(El, Ev, Esen, Eaz, Einc, Ed)}]
    else:
        modules = [{'label': 'internode', 'ntop': ntop, 'length': El,
                    'visible_length': Ev, 'senesced_length': Esen,
                    'diameter': Ed, 'azimuth': Eaz, 'inclination': Einc,
                    'elements': internode_elements(El, Ev, Esen, Eaz, Einc,
                                                   Ed)},
                   {'label': 'sheath', 'ntop': ntop, 'length': Gl,
                    'visible_length': Gv, 'senesced_length': Gsen,
                    'diameter': Gd, 'azimuth': Gaz, 'inclination': Ginc,
                    'elements': sheath_elements(Gl, Gv, Gsen, Gaz, Ginc, Gd)},
                   {'label': 'blade', 'ntop': ntop, 'length': Ll,
                    'rolled_length': Lr, 'visible_length': Lv,
                    'senesced_length': Lsen, 'n_sect': Lsect,
                    'shape_mature_length': L_shape, 'shape_max_width': Lw_shape,
                    'shape_key': shape_key, 'inclination': Linc,
                    'elements': blade_elements(Lsect, Ll, Lv, Lr, Lsen, L_shape,
                                               Lw_shape, shape_key, leaves,
                                               split=split)}]

        try:
            exposition = float(exposition)
            lifetime = float(lifetime)
            age = float(age)
            is_ligulated = int(is_ligulated)
            modules[2].update(
                {'exposition': exposition, 'lifetime': lifetime, 'age': age,
                 'is_ligulated': is_ligulated})
        except ValueError:
            pass

        if elongation:
            modules[0]['elongation_curve'] = {
                'x': [elongation['endleaf'], elongation['endE']], 'y': [0, El]}
            modules[1]['elongation_curve'] = {
                'x': [elongation['endBlade'], elongation['endleaf']],
                'y': [0, Gl]}
            modules[2]['elongation_curve'] = {
                'x': [elongation['startleaf'], elongation['endBlade']],
                'y': [0, Ll]}

    return modules


def get_component(components, index):
    component = components[index]
    elements = component['elements']
    properties = dict(component)
    del properties['elements']
    return properties, elements


def find_label(label, g, constrained_in=0):
    """
    Find the vertex ids identified by label in g if it exists.
    Args:
        label: (str) the label to be matched
        g: a MTG
        constrained_in: a vid defining the complex the label is search


    Returns: (list of int) the vertex id of matching labels

    """
    ci = constrained_in
    labels = g.property('label')
    if ci == 0:
        return [k for k in g.vertices_iter() if labels.get(k, '') == label]
    else:
        # Hack car iter_mtg peut recycler sur le reste du mtg
        cscale = g.scale(ci)
        return [k for k in iter_mtg(g, ci) if (labels.get(k, '') == label) and (
        g.complex_at_scale(k, cscale) == ci)]


def find_plants(g):
    """ the vid of all plants in g

    Args:
        g: an mtg

    Returns: (list of int) the vid of plants in g

    """
    vids = g.component_roots_at_scale(g.root, 1)
    if len(vids) > 0:
        labels = g.property('label')
        return [vid for vid in vids if labels[vid].startswith('plant')]
    else:
        return vids


def find_metamers(g, plant_number=1, axe_label='MS'):
    """ the vid of all metamers beared by axe indentified by its label

    Args:
        g: an mtg
        plant_number: (int) the number of the plant
        axe_label: (str) the label of the axe

    Returns:
        vid_plant (int) the vid of the plant bearing the metamers(int)
        vid_axe (int) the vid of the axe bearing the metamers
        vid_metamers (list of int) the vid of metamers

    """

    plant_label = 'plant' + str(plant_number)
    vid_plant = find_label(plant_label, g)
    if len(vid_plant) > 0:
        vid_plant = vid_plant[0]
    else:
        raise ValueError('plant ' + plant_label + ' not found in g')

    plant_axes = g.components(vid_plant)
    vid_axe = [vid for vid in find_label(axe_label, g) if vid in plant_axes]
    if len(vid_axe) > 0:
        vid_axe = vid_axe[0]
    else:
        raise ValueError(
            'axis ' + axe_label + ' not found in g on plant ' + plant_label)

    vids = g.components(vid_axe)
    if len(vids) > 0:
        labels = g.property('label')
        metamers = [vid for vid in vids if labels[vid].startswith('metamer')]
    else:
        metamers = []
    return vid_plant, vid_axe, metamers


def add_plant(g, plant_number=None, plant_properties=None, axis_properties=None,
              metamer_properties=None, collar_properties=None):
    """ Add a plant identified by its number to an mtg.
    The plant is created with 5 sub-scales: plant / mainstem /
    metamer0 / collar / baseElement
    If the plant already exists, nothing is added

    Args:
        g: an mtg representing the canopy
        plant_number: (int) the number of the plant to add.
        If None a plant_number will be max(existing_plant_numbers) + 1
        plant_properties: a dict of properties associated to the plant vertex
        axis_properties: a dict of properties associated to the mainstem
        metamer_properties: a dict of properties associated to metamer0
        collar_properties:  a dict of properties associated to collar

    Returns: the vid of the plant created or found

    """
    plant_properties = plant_properties or {}
    axis_properties = axis_properties or {}
    metamer_properties = metamer_properties or {}
    collar_properties = collar_properties or {}

    if plant_number is None:
        labels = g.property('label')
        plant_labels = [labels[vid] for vid in
                        g.component_roots_at_scale(g.root, scale=1)]
        if len(plant_labels) > 0:
            plant_number = max(
                map(lambda x: int(x.split('plant')[1]), plant_labels)) + 1
        else:
            plant_number = 1

    label = 'plant' + str(plant_number)
    found = find_label(label, g)
    if len(found) > 0:
        return found[0]
    else:
        vid_plant = g.add_component(g.root, label=label, edge_type='/',
                                    **plant_properties)
        vid_axe = g.add_component(vid_plant, edge_type='/', label='MS',
                                  **axis_properties)
        vid_metamer = g.add_component(vid_axe, edge_type='/', label='metamer0',
                                      **metamer_properties)
        vid_organ = g.add_component(vid_metamer, edge_type='/', label='collar',
                                    **collar_properties)
        vid_base = g.add_component(vid_organ, edge_type='/',
                                   label='baseElement')
        vid_top = g.add_component(vid_organ, edge_type='/', label='topElement')
        g.add_child(vid_base, child=vid_top, edge_type='<')
        return vid_plant


def add_vegetative_metamer(g, plant_number=1, axe_label='MS',
                           metamer_properties=None, internode_properties=None,
                           sheath_properties=None, blade_properties=None):
    """ Add a vegetative metatmer at the top of an axe

    Args:
        g: the MTG
        plant_number: (int) the plant number
        axe_label: the label of the axe bearing the new metamer
        metamer_properties: a dict of properties associated to the metamer
        internode_properties: a dict of properties associated to the internode
        sheath_properties: a dict of properties associated to the sheath
        blade_properties: a dict of properties associated to the blade

    Returns:
        the vid of the new metamer
    """

    metamer_properties = metamer_properties or {}
    internode_properties = internode_properties or {}
    sheath_properties = sheath_properties or {}
    blade_properties = blade_properties or {}

    vid_plants, vid_axe, metamers = find_metamers(g, plant_number, axe_label)
    num_metamer = len(metamers)  # start at zero
    label = 'metamer' + str(num_metamer)
    # to check with christophe : I want the 'tip'
    vid_parent_metamer = max(g.components(vid_axe))
    vid_metamer = g.add_component(vid_axe, edge_type='/', label=label,
                                  **metamer_properties)
    g.add_child(vid_parent_metamer, child=vid_metamer, edge_type='<')
    # add organs
    # parent = internode or collar = first organ of the preceding metamer
    vid_parent_organ = min(g.components(vid_parent_metamer))
    # Top of the parent organ
    vid_parent_elt = find_label('topElement', g, vid_parent_organ)[0]
    # add internode
    new_organ = g.add_component(vid_metamer, label='internode', edge_type='/',
                                **internode_properties)
    base_elt = g.add_component(new_organ, label='baseElement', edge_type='/')
    top_elt = g.add_component(new_organ, label='topElement', edge_type='/')
    vid_parent_organ = g.add_child(vid_parent_organ, child=new_organ,
                                   edge_type='<')
    g.add_child(vid_parent_elt, child=base_elt, edge_type='<')
    vid_parent_elt = g.add_child(base_elt, child=top_elt, edge_type='<')
    # sheath
    new_organ = g.add_component(vid_metamer, label='sheath', edge_type='/',
                                **sheath_properties)
    base_elt = g.add_component(new_organ, label='baseElement', edge_type='/')
    top_elt = g.add_component(new_organ, label='topElement', edge_type='/')
    vid_parent_organ = g.add_child(vid_parent_organ, child=new_organ,
                                   edge_type='+')
    g.add_child(vid_parent_elt, child=base_elt, edge_type='+')
    vid_parent_elt = g.add_child(base_elt, child=top_elt, edge_type='<')
    # blade
    new_organ = g.add_component(vid_metamer, label='blade', edge_type='/',
                                **blade_properties)
    base_elt = g.add_component(new_organ, label='baseElement', edge_type='/')
    top_elt = g.add_component(new_organ, label='topElement', edge_type='/')
    g.add_child(vid_parent_organ, child=new_organ, edge_type='<')
    g.add_child(vid_parent_elt, child=base_elt, edge_type='<')
    g.add_child(base_elt, child=top_elt, edge_type='<')

    return vid_metamer


def insert_elements(g, vid_organ, elements):
    """ Insert elements between base and top element of an organ

    Args:
        g: the mtg
        vid_organ: the vertex id of the targetd organ
        elements: a list of dict with element properties

    Returns:
        a list of vid of the inserted elements
    """

    inserted = []
    before = find_label('topElement', g, vid_organ)[0]
    for element in reversed(elements):
        before = g.insert_parent(before, edge_type='<', **element)
        inserted.append(before)
    return inserted


def add_axe(g, label, plant_number=1, axis_properties=None,
            metamer_properties=None, collar_properties=None):
    """ Add an axe identified by its label to a plant identifed by its number.
    The axe is created with 4 sub-scales axe /metamer0 / collar / baseElement
    If the axe already exists, nothing is added
    If bearing metamers are missing they are added

    Args:
        g: an mtg representing the canopy
        label: (str) a string identifying the axe (eg 'T1' for an axe beared by
        metamer 1 of main stem, 'T0.1' for the axe beared by metamer1 of T0)
        plant_number: (int) the number of the plant bering the axe.
        axis_properties: a dict of properties associated to the mainstem
        metamer_properties: a dict of properties associated to metamer0
        collar_properties:  a dict of properties associated to collar

    Returns: the vid of the axe created or found

    """
    plant_label = 'plant' + str(plant_number)
    vid_plant = find_label(plant_label, g)
    if len(vid_plant) > 0:
        vid_plant = vid_plant[0]
    else:
        raise ValueError('plant ' + plant_label + ' not found in g')

    axes = g.components(vid_plant)
    found = [vid for vid in find_label(label, g) if vid in axes]
    if len(found) > 0:
        return found[0]

    axis_properties = axis_properties or {}
    metamer_properties = metamer_properties or {}
    collar_properties = collar_properties or {}

    labels = g.property('label')
    axe_code = label.split('T')[1].split('.')

    if len(axe_code) > 1:
        parent_axe_label = 'T' + '.'.join(axe_code[:-1])
    else:
        parent_axe_label = 'MS'

    # recursively add parent axes if needed
    found = [vid for vid in find_label(parent_axe_label, g) if vid in axes]
    if len(found) == 0:
        add_axe(g, parent_axe_label, plant_number)

    vid_axe = g.add_component(vid_plant, edge_type='/', label=label,
                              **axis_properties)
    vid_metamer = g.add_component(vid_axe, edge_type='/', label='metamer0',
                                  **metamer_properties)
    vid_organ = g.add_component(vid_metamer, edge_type='/', label='collar',
                                **collar_properties)
    vid_base_elt = g.add_component(vid_organ, edge_type='/',
                                   label='baseElement')
    vid_top_elt = g.add_component(vid_organ, edge_type='/', label='topElement')

    vid_plant, vid_parent_axe, parent_metamers = find_metamers(g, plant_number,
                                                               parent_axe_label)
    bearing_metamer = int(axe_code[-1])
    last_metamer = len(parent_metamers) - 1
    if bearing_metamer > last_metamer:
        for i in range(last_metamer + 1, bearing_metamer + 1):
            new = add_vegetative_metamer(g, plant_number,
                                         labels[vid_parent_axe])
            parent_metamers.append(new)

    vid_parent_metamer = sorted(parent_metamers)[bearing_metamer]
    # parent organ = internode = first organ of the preceding metamer
    vid_parent_organ = min(g.components(vid_parent_metamer))
    # base of the parent organ = first element
    vid_parent_elt = find_label('baseElement', g, vid_parent_organ)[0]

    g.add_child(vid_parent_axe, child=vid_axe, edge_type='+')
    g.add_child(vid_parent_metamer, child=vid_metamer, edge_type='+')
    g.add_child(vid_parent_organ, child=vid_organ, edge_type='+')
    g.add_child(vid_parent_elt, child=vid_base_elt, edge_type='+')
    g.add_child(vid_base_elt, child=vid_top_elt, edge_type='<')

    return vid_axe


def update_organ_elements(g, leaves=None, split=False):
    """ Set / update organ elements

    Args:
        g: an adel mtg
        leaves: a leaf shape database
        split: (bool) flag trigering the separation between senescent and green
        part of an organ

    Returns:

    """
    labels = g.property('label')
    length = g.property('length')
    visible_length = g.property('visible_length')
    rolled_length = g.property('rolled_length')
    senesced_length = g.property('senesced_length')
    azimuth = g.property('azimuth')
    inclination = g.property('inclination')
    diameter = g.property('diameter')
    sectors = g.property('n_sect')
    shape_mature_length = g.property('shape_mature_length')
    shape_max_width = g.property('shape_max_width')
    shape_key = g.property('shape_key')

    for organ in g.vertices(scale=4):
        if labels[organ].startswith('internode'):
            elts = internode_elements(length[organ], visible_length[organ],
                                      senesced_length[organ], azimuth[organ],
                                      inclination[organ], diameter[organ],
                                      split=split)
        elif labels[organ].startswith('sheath'):
            elts = sheath_elements(length[organ], visible_length[organ],
                                   senesced_length[organ], azimuth[organ],
                                   inclination[organ], diameter[organ],
                                   split=split)
        elif labels[organ].startswith('blade'):
            elts = blade_elements(sectors[organ], length[organ],
                                  visible_length[organ], rolled_length[organ],
                                  senesced_length[organ],
                                  shape_mature_length[organ],
                                  shape_max_width[organ], shape_key[organ],
                                  leaves=leaves, split=split)
        else:
            elts = []

        if len(elts) > 0:
            if len(g.components(organ)) == 2:
                insert_elements(g, organ, elts)
            else:
                for elt in elts:
                    label = elt.pop('label')
                    vid_elt = find_label(label, g, organ)[0]
                    for k in elt:
                        g.property(k)[vid_elt] = elt[k]

    return g


def new_mtg_factory(parameters, metamer_factory=None, leaf_sectors=1,
                    leaves=None, stand=None, axis_dynamics=None,
                    add_elongation=False,
                    topology=('plant', 'axe_id', 'numphy'), split=False,
                    aborting_tiller_reduction=1.0, leaf_db=None):
    """A 'clone' of mtg_factory that uses mtg_edition functions
    """

    if leaf_db is not None:
        raise AdelDeprecationError(
            'leaf_db argument is deprecated, use leaves argument instead')

    if leaves is None:
        dynamic_leaf_db = False
    else:
        dynamic_leaf_db = leaves.dynamic

    g = MTG()

    # buffers
    # for detection of newplant/newaxe
    prev_plant = 0
    prev_axe = -1

    dp = parameters
    nrow = len(dp['plant'])

    for i in range(nrow):
        plant, num_metamer = [int(convert(dp.get(x)[i], undef=None)) for x in
                              [topology[e] for e in [0, 2]]]
        axe = dp.get(topology[1])[i]
        args = properties_from_dict(dp, i, exclude=topology)

        # Add plant if new
        if plant != prev_plant:
            if axe != 'MS':
                raise ValueError('Main stem is expected first when a new plant '
                                 'is declared')

            position, azimuth = (0, 0, 0), 0
            if stand and len(stand) >= plant:
                position, azimuth = stand[plant - 1]
            plant_properties = {'position': position, 'azimuth': azimuth,
                'refplant_id': args.get('refplant_id')}

            timetable = None
            if axis_dynamics:
                timetable = axis_dynamics[str(plant)][str(axe)]
            mainstem_properties = dict(timetable=timetable,
                                       HS_final=args.get('HS_final'),
                                       nff=args.get('nff'),
                                       hasEar=args.get('hasEar'),
                                       azimuth=args.get('az_insertion'))

            add_plant(g, plant_number=plant, plant_properties=plant_properties,
                      axis_properties=mainstem_properties)

        # Add axis
        if axe != prev_axe and axe != 'MS':
            timetable = None
            if axis_dynamics:
                timetable = axis_dynamics[str(plant)][str(axe)]
            axe_properties = dict(timetable=timetable,
                                  HS_final=args.get('HS_final'),
                                  nff=args.get('nff'),
                                  hasEar=args.get('hasEar'),
                                  azimuth=args.get('az_insertion'))
            add_axe(g, axe, plant, axis_properties=axe_properties)

        # Add metamer
        assert num_metamer > 0
        # args are added to metamers only if metamer_factory is none,
        # otherwise compute metamer components
        components = []
        if metamer_factory:
            xysr_key = None
            if leaves is not None and 'LcType' in args and 'LcIndex' in args:
                lctype = int(args['LcType'])
                lcindex = int(args['LcIndex'])
                if lctype != -999 and lcindex != -999:
                    age = None
                    if dynamic_leaf_db:
                        age = float(args[
                                        'rph']) - 0.3  # age_db = HS - rank + 1 = ph - 1.3 - rank +1 = rph - .3
                        if age != 'NA':
                            age = max(0, int(float(age)))
                    xysr_key = leaves.get_leaf_key(lctype, lcindex, age)

            elongation = None
            if add_elongation:
                startleaf = -.4
                endleaf = 1.6
                stemleaf = 1.2
                startE = endleaf
                endE = startE + (endleaf - startleaf) / stemleaf
                endBlade = endleaf
                if args['Gl'] > 0:
                    endBlade = args['Ll'] / args['Gl'] * (endleaf - startleaf)
                elongation = {'startleaf': startleaf, 'endBlade': endBlade,
                              'endleaf': endleaf, 'endE': endE}
            if not 'ntop' in args:
                args.update({'ntop': None})
            if not 'Gd' in args:
                args.update({'Gd': 0.19})
            args.update({'split': split})
            if args.get('HS_final') < args.get('nff'):
                for what in (
                        'Ll', 'Lv', 'Lr', 'Lsen', 'L_shape', 'Lw_shape', 'Gl',
                        'Gv', 'Gsen', 'Gd', 'El', 'Ev', 'Esen', 'Ed'):
                    args.update(
                        {what: args.get(what) * aborting_tiller_reduction})
            components = metamer_factory(Lsect=leaf_sectors, shape_key=xysr_key,
                                         elongation=elongation, leaves=leaves,
                                         **args)
            args = {'L_shape': args.get('L_shape')}
        #
        metamer_properties = args
        internode, sheath, blade = None, None, None
        elts = {k: [] for k in ('internode', 'sheath', 'blade')}
        if len(components) > 0:
            internode, sheath, blade = components
            internode.pop('label')
            sheath.pop('label')
            blade.pop('label')
            elts['internode'] = internode.pop('elements')
            elts['sheath'] = sheath.pop('elements')
            elts['blade'] = blade.pop('elements')

        vid_metamer = add_vegetative_metamer(g, plant, axe,
                                             metamer_properties=metamer_properties,
                                             internode_properties=internode,
                                             sheath_properties=sheath,
                                             blade_properties=blade)
        for organ in g.components(vid_metamer):
            label = g.property('label')[organ]
            if label in elts and len(elts[label]) > 0:
                insert_elements(g, organ, elts[label])
    return fat_mtg(g)


def mtg_factory(parameters, metamer_factory=None, leaf_sectors=1, leaves=None,
                stand=None, axis_dynamics=None, add_elongation=False,
                topology=['plant', 'axe_id', 'numphy'], split=False,
                aborting_tiller_reduction=1.0, leaf_db=None):
    """ Construct a MTG from a dictionary of parameters.

    The dictionary contains the parameters of all metamers in the stand (topology + properties).
    metamer_factory is a function that build metamer properties and metamer elements from parameters dict.
    leaf_sectors is an integer giving the number of LeafElements per Leaf blade
    leaves is an instance of adel.geometric_elements.Leaves class
    stand is a list of tuple (xy_position_tuple, azimuth) of plant positions
    axis_dynamics is a 3 levels dict describing axis dynamic. 1st key level is plant number, 2nd key level is axis number, and third ky level are labels of values (n, tip, ssi, disp)
    topology is the list of key names used in parameters dict for plant number, axe number and metamer number
    aborting_tiller_reduction is a scaling factor applied to reduce all dimensions of organs of tillers that will abort

    Axe number 0 is compulsory
  
    """

    if leaf_db is not None:
        raise AdelDeprecationError(
            'leaf_db argument is deprecated, use leaves argument instead')

    if leaves is None:
        dynamic_leaf_db = False
    else:
        dynamic_leaf_db = leaves.dynamic

    g = MTG()

    # buffers
    # for detection of newplant/newaxe
    prev_plant = 0
    prev_axe = -1
    # current vid
    vid_plant = -1
    vid_axe = -1
    vid_metamer = -1
    vid_node = -1
    vid_elt = -1
    # vid of top of stem nodes and elements
    vid_topstem_node = -1
    vid_topstem_element = -1
    # vid of plant main stem (axe0)
    vid_main_stem = -1
    # buffer for the vid of main stem anchor points for the first metamer, node and element of tillers
    metamers = []
    nodes = []
    elts = []

    dp = parameters
    nrow = len(dp['plant'])

    for i in range(nrow):
        plant, num_metamer = [int(convert(dp.get(x)[i], undef=None)) for x in
                              [topology[e] for e in [0, 2]]]
        axe = dp.get(topology[1])[i]
        mspos = int(convert(dp.get('ms_insertion')[i], undef=None))
        args = properties_from_dict(dp, i, exclude=topology)
        # Add plant if new
        if plant != prev_plant:
            label = 'plant' + str(plant)
            position = (0, 0, 0)
            azimuth = 0
            if stand and len(stand) >= plant:
                position, azimuth = stand[plant - 1]
            vid_plant = g.add_component(g.root, label=label, edge_type='/',
                                        position=position, azimuth=azimuth,
                                        refplant_id=args.get('refplant_id'))
            # reset buffers
            prev_axe = -1
            vid_axe = -1
            vid_metamer = -1
            vid_node = -1
            vid_elt = -1
            vid_topstem_node = -1
            vid_topstem_element = -1
            vid_main_stem = -1
            metamers = []
            nodes = []
            elts = []

        # Add axis
        if axe != prev_axe:
            label = ''.join(axe.split('.'))
            timetable = None
            if axis_dynamics:
                timetable = axis_dynamics[str(plant)][str(axe)]
            if axe == 'MS':
                vid_axe = g.add_component(vid_plant, edge_type='/', label=label,
                                          timetable=timetable,
                                          HS_final=args.get('HS_final'),
                                          nff=args.get('nff'),
                                          hasEar=args.get('hasEar'),
                                          azimuth=args.get('az_insertion'))
                vid_main_stem = vid_axe
            else:
                vid_axe = g.add_child(vid_main_stem, edge_type='+', label=label,
                                      timetable=timetable,
                                      HS_final=args.get('HS_final'),
                                      nff=args.get('nff'),
                                      hasEar=args.get('hasEar'),
                                      azimuth=args.get('az_insertion'))

        # Add metamer
        assert num_metamer > 0
        # args are added to metamers only if metamer_factory is none, otherwise compute metamer components
        components = []
        if metamer_factory:
            xysr_key = None
            if leaves is not None and 'LcType' in args and 'LcIndex' in args:
                lctype = int(args['LcType'])
                lcindex = int(args['LcIndex'])
                if lctype != -999 and lcindex != -999:
                    age = None
                    if dynamic_leaf_db:
                        age = float(args[
                                        'rph']) - 0.3  # age_db = HS - rank + 1 = ph - 1.3 - rank +1 = rph - .3
                        if age != 'NA':
                            age = max(0, int(float(age)))
                    xysr_key = leaves.get_leaf_key(lctype, lcindex, age)

            elongation = None
            if add_elongation:
                startleaf = -.4
                endleaf = 1.6
                stemleaf = 1.2
                startE = endleaf
                endE = startE + (endleaf - startleaf) / stemleaf
                endBlade = endleaf
                if args['Gl'] > 0:
                    endBlade = args['Ll'] / args['Gl'] * (endleaf - startleaf)
                elongation = {'startleaf': startleaf, 'endBlade': endBlade,
                              'endleaf': endleaf, 'endE': endE}
            if not 'ntop' in args:
                args.update({'ntop': None})
            if not 'Gd' in args:
                args.update({'Gd': 0.19})
            args.update({'split': split})
            if args.get('HS_final') < args.get('nff'):
                for what in (
                        'Ll', 'Lv', 'Lr', 'Lsen', 'L_shape', 'Lw_shape', 'Gl',
                        'Gv', 'Gsen', 'Gd', 'El', 'Ev', 'Esen', 'Ed'):
                    args.update(
                        {what: args.get(what) * aborting_tiller_reduction})
            components = metamer_factory(Lsect=leaf_sectors, shape_key=xysr_key,
                                         elongation=elongation, leaves=leaves,
                                         **args)
            args = {'L_shape': args.get('L_shape')}
        #
        label = 'metamer' + str(num_metamer)
        new_metamer = g.add_component(vid_axe, edge_type='/', label=label,
                                      **args)
        if axe == 'MS' and num_metamer == 1:
            vid_metamer = new_metamer
        elif num_metamer == 1:
            # add the edge with the bearing metamer on main stem
            vid_metamer = metamers[mspos - 1]
            vid_metamer = g.add_child(vid_metamer, child=new_metamer,
                                      edge_type='+')
        else:
            vid_metamer = g.add_child(vid_metamer, child=new_metamer,
                                      edge_type='<')

        # add metamer components, if any           
        if len(components) > 0:
            # deals with first component (internode) and first element 
            node, elements = get_component(components, 0)
            element = elements[0]
            new_node = g.add_component(vid_metamer, edge_type='/', **node)
            new_elt = g.add_component(new_node, edge_type='/', **element)
            if axe == 'MS' and num_metamer == 1:  # root of main stem
                vid_node = new_node
                vid_elt = new_elt
            elif num_metamer == 1:  # root of tiller
                vid_node = nodes[mspos - 1]
                vid_node = g.add_child(vid_node, child=new_node, edge_type='+')
                vid_elt = elts[mspos - 1]
                vid_elt = g.add_child(vid_elt, child=new_elt, edge_type='+')
            else:
                vid_node = g.add_child(vid_topstem_node, child=new_node,
                                       edge_type='<')
                vid_elt = g.add_child(vid_topstem_element, child=new_elt,
                                      edge_type='<')
            # add other elements of first component (the internode)
            for i in range(1, len(elements)):
                element = elements[i]
                vid_elt = g.add_child(vid_elt, edge_type='<', **element)
            vid_topstem_node = vid_node
            vid_topstem_element = vid_elt  # last element of internode

            # add other components   
            for i in range(1, len(components)):
                node, elements = get_component(components, i)
                if node['label'] == 'sheath':
                    edge_type = '+'
                else:
                    edge_type = '<'
                vid_node = g.add_child(vid_node, edge_type=edge_type, **node)
                element = elements[0]
                new_elt = g.add_component(vid_node, edge_type='/', **element)
                vid_elt = g.add_child(vid_elt, child=new_elt,
                                      edge_type=edge_type)
                for j in range(1, len(elements)):
                    element = elements[j]
                    vid_elt = g.add_child(vid_elt, edge_type='<', **element)

                    # update buffers
        if axe == 'MS':
            metamers.append(vid_metamer)
            if len(components) > 0:
                nodes.append(vid_topstem_node)
                elts.append(vid_topstem_element)
        prev_plant = plant
        prev_axe = axe

    return fat_mtg(g)


def update_elements(organ, leaves=None):
    if organ.label.startswith('blade'):
        elements = blade_elements(organ.n_sect, organ.length,
                                  organ.visible_length, organ.rolled_length,
                                  organ.senesced_length,
                                  organ.shape_mature_length,
                                  organ.shape_max_width, organ.shape_key,
                                  leaves=leaves)
        for i, e in enumerate(organ.components()):
            for k in elements[i]:
                exec "e.%s = elements[i]['%s']" % (k, k)


def update_plant(plant, time):
    """ update phenology of plant axes """
    plant.time = time


def update_axe(axe):
    """ update phenology on axes """
    if 'timetable' in axe.properties():
        axe.phyllochronic_time = numpy.interp(axe.complex().time,
                                              axe.timetable['tip'],
                                              axe.timetable['n'])
        # print 'axe %s phyllochronic time:%f'%(axe.label,axe.phyllochronic_time)


def update_organ(organ, h_whorl=0):
    rank = int(organ.complex().index())
    axe = organ.complex_at_scale(2)
    rph = axe.phyllochronic_time - rank
    length = organ.length
    vlength = organ.visible_length
    if 'elongation_curve' in organ.properties():
        organ.length = numpy.interp(rph, organ.elongation_curve['x'],
                                    organ.elongation_curve['y'])
    organ.visible_length = organ.length - h_whorl
    update_elements(organ)
    organ.dl = organ.length - length
    organ.dl_visible = organ.visible_length - vlength


def update_organ_from_table(organ, metamer, oldmetamer):
    neworg = metamer[organ.label]
    oldorg = oldmetamer[organ.label]
    new_elts = neworg.pop('elements')
    old_elts = oldorg.pop('elements')
    for k in neworg:
        if k is not 'shape_xysr':
            exec "organ.%s = neworg['%s']" % (k, k)
    for i, e in enumerate(organ.components()):
        has_area = False
        for k in new_elts[i]:
            if k in ['area', 'green_area', 'senesced_area']:
                exec "e.%s += (new_elts[i]['%s'] - old_elts[i]['%s'])" % (
                    k, k, k)
                has_area = True
            else:
                exec "e.%s = new_elts[i]['%s']" % (k, k)
        # control senescence (in case of acceleration by an other process)
        if has_area:
            if (e.green_area + e.senesced_area) > e.area:
                e.green_area = 0
                e.senesced_area = e.area


def mtg_update_at_time(g, time):
    """ Compute plant state at a given time according to dynamical parameters found in mtg
    """
    for pid in g.component_roots_at_scale_iter(g.root, scale=1):
        p = g.node(pid)
        update_plant(p, time)
        hw = {'0': 0}
        for a in p.components_at_scale(2):
            update_axe(a)
            numaxe = int(a.index())
            hwhorl = hw[str(numaxe)]
            for o in a.components_at_scale(4):
                update_organ(o, hwhorl)
                # update whorl
                if o.label.startswith('internode'):
                    if o.inclination < 0:
                        hwhorl = 0  # redressement
                    else:
                        hwhorl = max(0, hwhorl - o.length)
                elif o.label.statswith('sheath'):
                    hwhorl += o.visible_length
                    # memorise main stem values
                    if numaxe == 0:
                        hw[o.complex().index()] = o.length


def mtg_update_from_table(g, cantable, old_cantable):
    """ Compute plant state at a given time according to dynamical parameters found in mtg
    """

    df = pandas.DataFrame(cantable)
    old_df = pandas.DataFrame(old_cantable)
    for pid in g.component_roots_at_scale_iter(g.root, scale=1):
        p = g.node(pid)
        for ax in p.components_at_scale(2):
            for m in ax.components_at_scale(3):
                nump = int(''.join(list(p.label)[5:]))
                numphy = int(''.join(list(m.label)[7:]))
                dm = df[(df['plant'] == nump) & (df['axe_id'] == ax.label) & (
                    df['numphy'] == numphy)]
                old_dm = old_df[(old_df['plant'] == nump) & (
                    old_df['axe_id'] == ax.label) & (
                                    old_df['numphy'] == numphy)]
                # G. Garin 02/08: Addition of the following condition
                if (len(dm) > 0):
                    dmd = dict(
                        [(k, v[0]) for k, v in dm.to_dict('list').iteritems()])
                    old_dmd = dict([(k, v[0]) for k, v in
                                    old_dm.to_dict('list').iteritems()])
                    blade = m.components_at_scale(4)[2]
                    dmd['xysr_shape'] = blade.shape_xysr
                    dmd['Lsect'] = blade.n_sect
                    old_dmd['xysr_shape'] = blade.shape_xysr
                    old_dmd['Lsect'] = blade.n_sect
                    met = adel_metamer(**dmd)
                    newmetamer = dict([(mm['label'], mm) for mm in met])
                    old_met = adel_metamer(**old_dmd)
                    oldmetamer = dict([(mm['label'], mm) for mm in old_met])
                    for o in m.components_at_scale(4):
                        update_organ_from_table(o, newmetamer, oldmetamer)


def adel_label(g, vid):
    label = 'undef'
    if g.scale(vid) == 5:
        organ = g.complex(vid)
        metamer = g.complex(organ)
        axe = g.complex(metamer)
        plant = g.complex(axe)
        label = '_'.join(
            [g.label(plant), g.label(axe), g.label(metamer), g.label(organ),
             g.label(vid)])
    return label


def adel_labels(g, scale=5):
    """ return a dict vid:adel_id
    """
    return {vid: adel_label(g, vid) for vid in g.vertices_iter(scale=scale)}


def adel_ids(g, scale=5):
    """ return a dict adel_id:vid
    """
    return {adel_label(g, vid): vid for vid in g.vertices_iter(scale=scale)}


def mtg_update(newg, g, refg):
    """ update newg with specific properties found in g only and update by increment compared to refg area-like properties
    """
    specific = set(g.property_names()) - set(newg.property_names())
    ids = adel_ids(g, scale=5)
    newids = adel_ids(newg, scale=5)
    common_labs = set(ids) & set(newids)

    for prop in specific:
        newg.add_property(prop)
        newprop = {newids[lab]: g.property(prop)[ids[lab]] for lab in
                   common_labs if ids[lab] in g.property(prop)}
        newg.property(prop).update(newprop)
        # g.remove_property(prop)#helps beeing compatible with ctypes objects

    for lab in common_labs:

        vid = ids[lab]
        newvid = newids[lab]
        if vid in g.property('area') and newvid in newg.property(
                'area') and vid in refg.property('area'):
            if newg.property('length')[newvid] > 0:
                dlength = max([0, newg.property('length')[newvid] -
                               refg.property('length')[vid]])
                darea = 0
                if dlength > 0:  # avoid changing area when length is stabilised
                    darea = max([0, newg.property('area')[newvid] -
                                 refg.property('area')[
                                     vid]])  # do not take into account negative variations du to rolling
                newarea = g.property('area')[vid] + darea

                # correction if last area estimation of area_sen is different from the one of area
                if newg.property('senesced_length')[newvid] >= \
                        newg.property('length')[newvid]:
                    newsen = newarea
                    newgreen = 0
                else:
                    dlength = max([0, newg.property('senesced_length')[newvid] -
                                   refg.property('senesced_length')[vid]])
                    dsen = 0
                    if dlength > 0:
                        dsen = max([0, newg.property('senesced_area')[newvid] -
                                    refg.property('senesced_area')[vid]])

                    newsen = min(
                        [newarea, g.property('senesced_area')[vid] + dsen])
                    newgreen = max([0, min([newarea - newsen,
                                            g.property('green_area')[vid] + (
                                                darea - dsen)])])

                newg.property('area')[newvid] = newarea
                newg.property('green_area')[newvid] = newgreen
                newg.property('senesced_area')[newvid] = newsen
            else:  # node is no longer there, specific property are removed
                for prop in specific:
                    if newvid in newg.property(prop):
                        newg.property(prop).pop(newvid)
    return newg


def move_properties(g_source, g_dest, filter_length=True, cleanup_source=True):
    """ Move properties present in g_source and not in g_dest into g_dest.
        if filter_length is True (default), properties attached to node whose length is zero are not transfered
    """
    specific = set(g_source.property_names()) - set(g_dest.property_names())
    ids = adel_ids(g_source, scale=5)
    newids = adel_ids(g_dest, scale=5)
    if filter_length:
        length = g_dest.property('length')
        newids = {lab: vid for lab, vid in newids.iteritems() if
                  length[vid] > 0}
    common_labs = set(ids) & set(newids)

    for prop in specific:
        g_dest.add_property(prop)
        newprop = {newids[lab]: g_source.property(prop)[ids[lab]] for lab in
                   common_labs if ids[lab] in g_source.property(prop)}
        g_dest.property(prop).update(newprop)
        if cleanup_source:
            g_source.remove_property(prop)


def exposed_areas(g):
    """ returns a Dataframe with all exposed (visible) areas of elements in g """
    data = {}
    what = ('length', 'area', 'green_length', 'green_area', 'senesced_length',
            'senesced_area')
    for vid in g.vertices_iter(scale=g.max_scale()):
        n = g.node(vid)
        if n.length > 0 and not n.label.startswith('Hidden'):
            numphy = int(''.join(list(n.complex().complex().label)[7:]))
            nf = n.complex().complex().complex().nff
            node_data = {
                'plant': n.complex().complex().complex().complex().label,
                'axe': n.complex().complex().complex().label, 'metamer': numphy,
                'organ': n.complex().label, 'vid': vid, 'ntop': nf - numphy + 1,
                'element': n.label,
                'refplant_id': n.complex().complex().complex().complex().refplant_id,
                'nff': nf, 'HS_final': n.complex().complex().complex().HS_final,
                'L_shape': n.complex().complex().L_shape}
            properties = n.properties()
            node_data.update({k: properties[k] for k in what})
            data[vid] = node_data
    df = pandas.DataFrame(data).T
    # hack
    df['d_basecol'] = 0
    return df


def exposed_areas2canS(exposed_areas):
    """ adaptor to convert new adel output to old adel output (canS-like) dataframe """
    d = exposed_areas
    if len(d) > 0:
        grouped = d.groupby(['plant', 'axe', 'metamer'], group_keys=False)

        def _metamer(sub):
            met = {'plant': sub.plant.values[0],
                   'refplant_id': sub.refplant_id.values[0],
                   'axe_id': sub.axe.values[0], 'nff': sub.nff.values[0],
                   'HS_final': sub.HS_final.values[0],
                   'numphy': sub.metamer.values[0], 'ntop': sub.ntop.values[0],
                   'L_shape': sub.L_shape.values[0],
                   'Lv': sub[sub.organ == 'blade'].length.sum(),
                   'Lvgreen': sub[sub.organ == 'blade'].green_length.sum(),
                   'Lvsen': sub[sub.organ == 'blade'].senesced_length.sum(),
                   'Slv': sub[sub.organ == 'blade'].area.sum(),
                   'Slvgreen': sub[sub.organ == 'blade'].green_area.sum(),
                   'Slvsen': sub[sub.organ == 'blade'].senesced_area.sum(),
                   'Gv': sub[sub.organ == 'sheath'].length.sum(),
                   'Gvgreen': sub[sub.organ == 'sheath'].green_length.sum(),
                   'Gvsen': sub[sub.organ == 'sheath'].senesced_length.sum(),
                   'SGv': sub[sub.organ == 'sheath'].area.sum(),
                   'SGvgreen': sub[sub.organ == 'sheath'].green_area.sum(),
                   'SGvsen': sub[sub.organ == 'sheath'].senesced_area.sum(),
                   'Ev': sub[sub.organ == 'internode'].length.sum(),
                   'Evgreen': sub[sub.organ == 'internode'].green_length.sum(),
                   'Evsen': sub[sub.organ == 'internode'].senesced_length.sum(),
                   'SEv': sub[sub.organ == 'internode'].area.sum(),
                   'SEvgreen': sub[sub.organ == 'internode'].green_area.sum(),
                   'SEvsen': sub[sub.organ == 'internode'].senesced_area.sum()}
            return pandas.DataFrame(met, index=[sub.index[0]])

        d = grouped.apply(_metamer)
        # hack
        d['d_basecol'] = 0
    return d

# to do

# varaibles sur element : area, senesced_area, senescenece_position (0,1 sur la feuille)
# generer dans alep.wheat.py un tableau axisdynamic
# la fonction grow
#
