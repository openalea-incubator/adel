"""
A place to develop candidates methods for mtg.py / topological builder

"""

#temporary import
from alinea.adel.mtg import convert,properties_from_dict

# import csv

from openalea.mtg import MTG, fat_mtg

import numpy as np
import pandas as pd



# try:
    # from openalea.mtg.traversal import *
# except:
    # pass

#from openalea.mtg.io import read_lsystem_string
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

def internode_elements(l,lvis, lsen, az, inc, split = False):
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
    try:
        lhide = max(l - lvis, 0.)
        lgreen = lvis - min(lsen,lvis)
        lsen = lvis - lgreen
        #
        is_green  = lgreen >= lsen
    except TypeError:
        pass
    if split :
        green_elt = {'label': 'StemElement', 'offset': lhide, 'length': lgreen, 'is_green': True, 'azimuth': az, 'inclination': inc}
        sen_elt = {'label': 'StemElement', 'offset': 0, 'length': lsen, 'is_green': False, 'azimuth': 0, 'inclination': 0}
        return [green_elt, sen_elt]
    else : 
        elt = {'label': 'StemElement', 'offset': lhide, 'length': lvis, 'l_sen' : lsen, 'is_green': is_green, 'azimuth': az, 'inclination': inc}
        return [elt]
    
def sheath_elements(l, lvis, lsen, az, inc, split = False):
    """ returns parameters of sheath elements (visible parts of the sheath).
    l is the length of the sheath
    lv is the visible length (senesced + green)
    lsen is the senescent apical length
    Fisrt elements is for green visible part of the sheath
    Second elements is for senesced visible part of the sheath
    """
    # same logic as internodes
    return internode_elements(l,lvis,lsen, az, inc, split)


def blade_elt_area(leaf, Lshape, Lwshape, sr_base, sr_top):
    """ surface of a blade element, positioned with two relative curvilinear absisca"""
    from scipy.integrate import simps
    from numpy import interp
    #try:
    #raise Exception("")
    S=0
    if leaf is not None:
        x,y,s,r = leaf
        sre = [sr for sr in zip(s,r) if (sr[0] > sr_base) & (sr[0] < sr_top)]
        # Temp G.Garin: 02/08/2013
        if len(sre)>0:
            se,re = zip(*sre)
            snew = [sr_base] + list(se) + [sr_top]
            rnew = [interp(sr_base,s,r)] + list(re) + [interp(sr_top,s,r)]
            S = simps(rnew,snew) * Lshape * Lwshape
        #print "S",S
    #except:
        #S = 0
    return S
    
def blade_elements(sectors, l, lvis, lrolled, lsen, Lshape, Lwshape, xysr_shape, split = False):
    """ return parameters of blade elements (visible parts of the blade).
    sectors is the number of sectors dividing pattern blade shape
    l is the length of the blade
    lvis is the visible length (senesced + green, rolled + flat)
    lrolled is the visible rolled length of the blade
    lsen is the senescent apical length
    Lshape is length of the blade used as a pattern shape
"""
    lhide = None
    lgreen = None
    lflat = None
    # s (on mature shape) at which leaf becomes flat and visible
    s_limvis = 1.
    #s (on mature shape) at which leaf becomes senescent
    s_limsen = 1.
    try:
        lhide = max(l - lvis, 0.)
        lflat = lvis - min(lrolled,lvis)
        lgreen = lflat - min(lsen,lflat)
        lsen = lflat - lgreen        
        s_limvis = Lshape - lflat
        s_limsen = Lshape - lsen
        #print(lgreen,lsen,s_limvis,s_limsen)
    except TypeError:
        pass

    # hidden part + rolled
    hidden_elt = {'label': 'StemElement', 'offset': lhide, 'length': lrolled, 'is_green': True}
    elements = [hidden_elt]
    ds = 0
    if Lshape:
        ds = float(Lshape) / sectors
    st = ds    
    for isect in range(sectors):
        ls_vis= 0
        ls_green = 0
        ls_sen = 0
        S_green = 0
        S_sen = 0
        position_senescence = 1
        srb_green = None
        srt_green = None
        srb_sen = None
        srt_sen = None
        try:
            ls_vis = max(0., st - s_limvis)
            if ls_vis > 0:
                sb_green = st - min(ds, ls_vis)
                st_green = min(st, max(sb_green,s_limsen))
                sb_sen = st_green
                st_sen= st
                ls_green = st_green - sb_green
                ls_sen = st_sen - sb_sen
                #print(sb_green,st_green,st_sen)
                if lflat > 0:
                    srb_green = float(sb_green - s_limvis) / lflat
                    srt_green = float(st_green - s_limvis) / lflat
                    srb_sen = float(sb_sen - s_limvis) / lflat
                    srt_sen = float(st_sen - s_limvis) / lflat
                    #print(srb_green,srt_green,srb_sen,srt_sen)
                    if ls_green > 0:
                        S_green = blade_elt_area(xysr_shape, Lshape, Lwshape, srb_green, srt_green)
                    if ls_sen > 0:
                        S_sen = blade_elt_area(xysr_shape, Lshape, Lwshape, srb_sen, srt_sen)
                    # attention a garder une position constante quand on utlise une feuille stresse
                    position_senescence = 1 - float(ls_sen) / Lshape
        except TypeError:# input is None
        #    print "passing"
            pass
        green_elt = {'label': 'LeafElement', 'length': ls_green, 'area': S_green, 'is_green': True,
                'srb': srb_green, 'srt': srt_green}
        sen_elt = {'label': 'LeafElement', 'length': ls_sen,'area': S_sen, 'is_green': False, 
                'srb': srb_sen, 'srt': srt_sen}
        elt = {'label': 'LeafElement', 'length': ls_sen + ls_green,'area': S_green + S_sen, 'green_area' : S_green, 'is_green': (ls_green > ls_sen), 
                'srb': srb_green, 'srt': srt_sen, 'position_senescence':position_senescence} 
        if split: 
            elements.extend([green_elt,sen_elt])
        else:
            elements.extend([elt])
        st += ds
    return elements
    
def adel_metamer(Ll=None, Lv=None, Lr=None, Lsen=None, L_shape=None, Lw_shape=None, xysr_shape=None, Linc=None, Laz=None, Lsect=1, Gl=None, Gv=None, Gsen=None, Gd=None, Ginc=None, El=None, Ev=None, Esen=None, Ed=None, Einc=None, elongation=None, **kwargs):
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
       - Gd : apparent diameter of the sheath
       - Ginc : relative inclination of the sheath
       - El: length of the internode (hidden + visible)
       - Ev: emerged length of the internode 
       - Esen: senescent length of the internode (hidden + visible)
       - Ed: diameter of the internode
       - Einc : relative inclination of the internode
 
    """
       #to do add diameter and Lrolled to blade
    Eaz = Laz
    Gaz = 0
    modules = [
        {'label': 'internode',
        'length': El,
        'visible_length': Ev,
        'senesced_length': Esen,
        'diameter' : Ed,
        'azimuth': Eaz,
        'inclination' : Einc,
        'elements' : internode_elements(El, Ev, Esen, Eaz, Einc)}, 
        {'label': 'sheath',
        'length': Gl,
        'visible_length': Gv,
        'senesced_length': Gsen,
        'diameter' : Gd,
        'azimuth' : Gaz,   
        'inclination' : Ginc,
        'elements': sheath_elements(Gl, Gv, Gsen, Gaz, Ginc)}, 
        {'label': 'blade',
        'length': Ll,
        'rolled_length': Lr,
        'diameter': Gd,
        'visible_length': Lv,
        'senesced_length': Lsen,
        'n_sect': Lsect,
        'shape_mature_length': L_shape,
        'shape_max_width' : Lw_shape,
        'shape_xysr': xysr_shape,
        'inclination' : Linc,
        'elements': blade_elements(Lsect, Ll, Lv, Lr, Lsen, L_shape, Lw_shape, xysr_shape)} 
    ]
    
    if elongation:
        
        modules[0]['elongation_curve'] = {'x': [elongation['endleaf'], elongation['endE']], 'y' : [0,El]}
        modules[1]['elongation_curve'] = {'x': [elongation['endBlade'], elongation['endleaf']], 'y' : [0,Gl]}
        modules[2]['elongation_curve'] = {'x': [elongation['startleaf'], elongation['endBlade']], 'y' : [0,Ll]}
    
    return modules
    
def get_component(components, index):
    component = components[index]
    elements = component['elements']
    properties = dict(component)
    del properties['elements']
    return properties, elements

    
def mtg_factory(parameters, metamer_factory=None, leaf_sectors=1, leaf_db = None, stand = None, axis_dynamics = None, add_elongation = False, topology = ['plant','axe_id','numphy']):
    """ Construct a MTG from a dictionary of parameters.

    The dictionary contains the parameters of all metamers in the stand (topology + properties).
    metamer_factory is a function that build metamer properties and metamer elements from parameters dict.
    leaf_sectors is an integer giving the number of LeafElements per Leaf blade
    leaf_db is xysr leaf shape database
    stand is a list of tuple (xy_position_tuple, azimuth) of plant positions
    axis_dynamics is a 3 levels dict describing axis dynamic. 1st key level is plant number, 2nd key level is axis number, and third ky level are labels of values (n, tip, ssi, disp)
    topology is the list of key names used in parameters dict for plant number, axe number and metamer number

    Axe number 0 is compulsory
  
    """
    
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
    #vid of plant main stem (axe0)
    vid_main_stem = -1
    # buffer for the vid of main stem anchor points for the first metamer, node and element of tillers
    metamers = []
    nodes = []
    elts = []

    dp = parameters
    nrow = len(dp['plant'])
    
    for i in range(nrow):
        plant, num_metamer = [int(convert(dp.get(x)[i],undef=None)) for x in [topology[e] for e in [0,2]]]        
        axe = dp.get(topology[1])[i]
        mspos = int(convert(dp.get('ms_insertion')[i],undef=None))
        # Add plant if new
        if plant != prev_plant:
            label = 'plant' + str(plant)
            position = (0,0,0)
            azimuth = 0
            if stand and len(stand) >= plant:
                position,azimuth = stand[plant-1]                
            vid_plant = g.add_component(g.root, label=label, edge_type='/', position = position, azimuth = azimuth)
            #reset buffers
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
            label = axe
            timetable = None
            if axis_dynamics:
                timetable = axis_dynamics[str(plant)][str(axe)]
            if axe == 'MS':
                vid_axe = g.add_component(vid_plant,edge_type='/',label=label, timetable=timetable)
                vid_main_stem = vid_axe
            else:
                vid_axe = g.add_child(vid_main_stem, edge_type='+',label=label, timetable=timetable)

        # Add metamer
        assert num_metamer > 0
        # args are added to metamers only if metamer_factory is none, otherwise compute metamer components
        args = properties_from_dict(dp,i,exclude=topology)
        components = []
        if metamer_factory:
            if leaf_db is not None:
                try:
                    xysr = leaf_db[str(int(args['LcType']))][int(args['LcIndex']) - 1]#R index starts at 1
                except KeyError:
                    xysr=leaf_db[leaf_db.keys()[0]][0]
            else:
                xysr = None
            
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
                elongation = {'startleaf' : startleaf , 'endBlade':endBlade, 'endleaf': endleaf, 'endE': endE}
            components = metamer_factory(Lsect = leaf_sectors, xysr_shape = xysr, elongation = elongation, **args)
            args={}
        #
        label = 'metamer'+str(num_metamer)
        new_metamer = g.add_component(vid_axe, edge_type='/', label = label, **args)
        if axe=='MS' and num_metamer==1:
            vid_metamer = new_metamer
        elif num_metamer == 1:
            # add the edge with the bearing metamer on main stem
            vid_metamer = metamers[mspos-1]
            vid_metamer =  g.add_child(vid_metamer, child=new_metamer, edge_type='+')
        else:
            vid_metamer = g.add_child(vid_metamer, child=new_metamer, edge_type='<')

        # add metamer components, if any           
        if len(components) > 0:
            # deals with first component (internode) and first element 
            node, elements = get_component(components,0)
            element = elements[0]
            new_node = g.add_component(vid_metamer, edge_type='/', **node)
            new_elt = g.add_component(new_node, edge_type='/', **element)
            if axe=='MS' and num_metamer==1: #root of main stem
                vid_node = new_node
                vid_elt =  new_elt                   
            elif num_metamer == 1: # root of tiller                   
                vid_node = nodes[mspos - 1]
                vid_node = g.add_child(vid_node, child = new_node, edge_type='+')
                vid_elt = elts[mspos - 1]
                vid_elt = g.add_child(vid_elt, child = new_elt, edge_type='+')
            else:
                vid_node = g.add_child(vid_topstem_node, child=new_node, edge_type='<')
                vid_elt = g.add_child(vid_topstem_element, child=new_elt, edge_type='<')
            #add other elements of first component (the internode)
            for i in range(1,len(elements)):
                element = elements[i]
                vid_elt = g.add_child(vid_elt, edge_type='<',**element)
            vid_topstem_node = vid_node
            vid_topstem_element = vid_elt #last element of internode  
            
            # add other components   
            for i in range(1,len(components)):
                node, elements = get_component(components,i)
                if node['label'] == 'sheath':
                    edge_type = '+'
                else:
                    edge_type = '<'
                vid_node = g.add_child(vid_node, edge_type=edge_type, **node)      
                element = elements[0]
                new_elt = g.add_component(vid_node, edge_type='/', **element)
                vid_elt = g.add_child(vid_elt, child=new_elt, edge_type=edge_type)
                for j in range(1,len(elements)):
                    element = elements[j]
                    vid_elt = g.add_child(vid_elt, edge_type='<',**element) 
                
        #update buffers 
        if axe == 'MS' :
            metamers.append(vid_metamer)
            if len(components) > 0:
                nodes.append(vid_topstem_node)
                elts.append(vid_topstem_element)
        prev_plant = plant
        prev_axe = axe
    
    return fat_mtg(g)


def update_elements(organ):
    if organ.label.startswith('blade'):
        rolled_length = 0
        elements =  blade_elements(organ.n_sect, organ.length, rolled_length, organ.visible_length, organ.senesced_length, organ.shape_mature_length)
        for i,e in enumerate(organ.components()):
            for k in elements[i]:
                exec "e.%s = elements[i]['%s']"%(k,k)

def update_plant(plant, time):
    """ update phenology of plant axes """
    plant.time = time
    
def update_axe(axe):
    """ update phenology on axes """
    if 'timetable' in axe.properties():
        axe.phyllochronic_time = np.interp(axe.complex().time, axe.timetable['tip'], axe.timetable['n'])
        #print 'axe %s phyllochronic time:%f'%(axe.label,axe.phyllochronic_time)
        
def update_organ(organ,h_whorl=0):
    rank = int(organ.complex().index())
    axe = organ.complex_at_scale(2)
    rph = axe.phyllochronic_time - rank
    length = organ.length
    vlength = organ.visible_length
    if 'elongation_curve' in organ.properties():
        organ.length = np.interp(rph, organ.elongation_curve['x'], organ.elongation_curve['y'])
    organ.visible_length = organ.length - h_whorl
    update_elements(organ)
    organ.dl = organ.length - length
    organ.dl_visible = organ.visible_length - vlength

def update_organ_from_table(organ,metamer, oldmetamer):
    neworg = metamer[organ.label]
    oldorg = oldmetamer[organ.label]
    new_elts = neworg.pop('elements')
    old_elts = oldorg.pop('elements')
    for k in neworg:
        if k is not 'shape_xysr':
            exec "organ.%s = neworg['%s']"%(k,k)
    for i,e in enumerate(organ.components()):
        for k in new_elts[i]:
            if k in ['area', 'green_area']:
                exec "e.%s += (new_elts[i]['%s'] - old_elts[i]['%s'])"%(k,k,k)
            else:
                exec "e.%s = new_elts[i]['%s']"%(k,k)

    

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
                update_organ(o,hwhorl)
                #update whorl
                if o.label.startswith('internode'):
                    if o.inclination < 0:
                        hwhorl = 0 #redressement
                    else:
                        hwhorl = max(0,hwhorl - o.length)
                elif o.label.statswith('sheath'):
                    hwhorl += o.visible_length
                    # memorise main stem values
                    if numaxe == 0:
                        hw[o.complex().index()] = o.length


def mtg_update_from_table(g, cantable, old_cantable): 
    """ Compute plant state at a given time according to dynamical parameters found in mtg
    """ 

    df = pd.DataFrame(cantable)
    old_df = pd.DataFrame(old_cantable)
    for pid in g.component_roots_at_scale_iter(g.root, scale=1):
        p = g.node(pid)
        for ax in p.components_at_scale(2):
            for m in ax.components_at_scale(3):
                nump = int(''.join(list(p.label)[5:]))
                numphy = int(''.join(list(m.label)[7:]))
                dm = df[(df['plant'] == nump) & (df['axe_id'] == ax.label) & (df['numphy'] == numphy)]
                old_dm = old_df[(old_df['plant'] == nump) & (old_df['axe_id'] == ax.label) & (old_df['numphy'] == numphy)]
                # G. Garin 02/08: Addition of the following condition
                if (len(dm) > 0):
                    dmd = dict([(k,v[0]) for k,v in dm.to_dict('list').iteritems()])
                    old_dmd = dict([(k,v[0]) for k,v in old_dm.to_dict('list').iteritems()])
                    blade = m.components_at_scale(4)[2]
                    dmd['xysr_shape'] = blade.shape_xysr
                    old_dmd['xysr_shape'] = blade.shape_xysr
                    met = adel_metamer(**dmd)
                    newmetamer = dict([(mm['label'],mm) for mm in met])
                    old_met = adel_metamer(**old_dmd)
                    oldmetamer = dict([(mm['label'],mm) for mm in old_met])
                    for o in m.components_at_scale(4):
                        update_organ_from_table(o,newmetamer,oldmetamer)

# to do

# varaibles sur element : area, senesced_area, senescenece_position (0,1 sur la feuille)
# generer dans alep.wheat.py un tableau axisdynamic
# la fonction grow
# 
