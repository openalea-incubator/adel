"""
A place to develop candidates methods for mtg.py

"""

#temporary import
from alinea.adel.mtg import convert,properties_from_dict

# import csv

from openalea.mtg import MTG, fat_mtg

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

def internode_elements(l,lvis, lsen):
    """ returns parameters of internode elements (visible parts of the internode).
    l is the length of the internode
    lv is the visible length (senesced + green)
    lsen is the senescent apical length
    Fisrt element is for the green visible part of the internode
    Second element is for senesced visible part of the internode
    """
    lhide = None
    lgreen = None
    try:
        lhide = max(l - lvis, 0.)
        lgreen = lvis - min(lsen,lvis)
        lsen = lvis - lgreen
    except TypeError:
        pass
    green_elt = {'label': 'StemElement', 'offset': lhide, 'length': lgreen, 'is_green': True}
    sen_elt = {'label': 'StemElement', 'offset': 0, 'length': lsen, 'is_green': False}
    return [green_elt, sen_elt]
    
def sheath_elements(l, lvis, lsen):
    """ returns parameters of sheath elements (visible parts of the sheath).
    l is the length of the sheath
    lv is the visible length (senesced + green)
    lsen is the senescent apical length
    Fisrt elements is for green visible part of the sheath
    Second elements is for senesced visible part of the sheath
    """
    # same logic as internodes
    return internode_elements(l,lvis,lsen)
    
def blade_elements(sectors, l, lvis, lsen, Lshape):
    """ return parameters of blade elements (visible parts of the blade).
    sectors is the number of sectors dividing pattern blade shape
    l is the length of the blade
    lvis is the visible length (senesced + green)
    lsen is the senescent apical length
    Lshape is length of the blade used as pattern shape
"""
    lhide = None
    lgreen = None
    #relative s (on mature shape) at which leaf becomes visible
    s_limvis = 1.
    #relative s (on mature shape) at which leaf becomes senescent
    s_limsen = 1.
    try:
        lhide = max(l - lvis, 0.)
        lgreen = lvis - min(lsen,lvis)
        lsen = lvis - lgreen        
        s_limvis = Lshape - lvis
        s_limsen = Lshape - lsen
        print(lgreen,lsen,s_limvis,s_limsen)
    except TypeError:
        pass

    # hidden part
    hidden_elt = {'label': 'StemElement', 'offset': lhide, 'length': 0, 'is_green': True}
    elements=[hidden_elt]
    ds = 0
    if Lshape:
        ds = float(Lshape) / sectors
    st = ds    
    for isect in range(sectors):
        ls_vis= 0
        ls_green = 0
        ls_sen = 0
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
                print(sb_green,st_green,st_sen)
                if lvis > 0:
                    srb_green = (sb_green - s_limvis) / lvis
                    srt_green = (st_green - s_limvis) / lvis
                    srb_sen = (sb_sen - s_limvis) / lvis
                    srt_sen = (st_sen - s_limvis) / lvis
        except TypeError:
            pass
        green_elt = {'label': 'LeafElement', 'length': ls_green, 'is_green': True,
                'srb': srb_green, 'srt': srt_green}
        sen_elt = {'label': 'LeafElement', 'length': ls_sen,'is_green': False, 
                'srb': srb_sen, 'srt': srt_sen} 
        elements.extend([green_elt,sen_elt])
        st += ds
    return elements
    
def adel_metamer(Ll=None, Lv=None, Lsen=None, L_shape=None, Lw_shape=None, xysr_shape=None, Linc=None, Laz=None, Lsect=1, Gl=None, Gv=None, Gsen=None, Gd=None, Ginc=None, El=None, Ev=None, Esen=None, Ed=None, Einc=None, *args):
    """ Contructs metamer elements for adel from parameters describing a static state.
    Parameters are : 
       - Ll : length of the blade
       - Lv : visible (emerged) length of blade (green + senesced)
       - Lsen : length of the senescent part of the blade (hidden + visible)       
       - L_shape : Mature length of the blade used to compute blade shape
       - Lw_shape : Maximal width of the blade used to compute blade shape
       - xysr_shape : a (x,y,s,r) tuple describing blade geometry
       - Linc : relative inclination of the base of leaf blade at the top of leaf sheath (deg) 
       - Laz : relative ? azimuth of the leaf
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
       
    modules = [
        {'label': 'internode',
        'length': El,
        'visible_length': Ev,
        'senesced_length': Esen,
        'diameter' : Ed, 
        'relative_inclination' : Einc,
        'elements' : internode_elements(El, Ev, Esen)}, 
        {'label': 'sheath',
        'length': Gl,
        'visible_length': Gv,
        'senesced_length': Gsen,
        'diameter' : Gd,
        'relative_inclination' : Ginc,
        'azimuth' : Laz,
        'elements': sheath_elements(Gl, Gv, Gsen)}, 
        {'label': 'blade',
        'length': Ll,
        'visible_length': Lv,
        'senesced_length': Lsen,
        'n_sect': Lsect,
        'shape_mature_length': L_shape,
        'shape_max_width' : Lw_shape,
        'shape_xysr': xysr_shape,
        'inclination' : Linc,
        'elements': blade_elements(Lsect, Ll, Lv, Lsen, L_shape)} 
    ]
    
    return modules
    
def get_component(components, index):
    component = components[index]
    elements = component['elements']
    properties = dict(component)
    del properties['elements']
    return properties, elements

    
def mtg_factory(parameters, metamer_factory=None, leaf_sectors=1, leaf_db = None, topology = ['plant','axe','numphy']):
    """ Construct a MTG from a dictionary of parameters.

    The dictionary contains the parameters of all metamers in the stand (topology + properties).
    metamer_factory is a function that build metamer properties and metamer elements from parameters dict.
    Sector is an integer giving the number of LeafElements per Leaf blade
    topology is the list of key names used in parameters dictfor plant number, axe numer and metamer number

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
    # buffer for the vid of main stem anchor points for the first metamer, node and element of tillers
    metamers = []
    nodes = []
    elts = []

    dp = parameters
    nrow = len(dp['plant'])
    
    for i in range(nrow):
        plant, axe, num_metamer = [int(convert(dp.get(x)[i],undef=None)) for x in topology]        

        # Add plant if new
        if plant != prev_plant:
            label = 'plant' + str(plant)
            vid_plant = g.add_component(g.root, label=label, edge_type='/')
            #reset buffers
            prev_axe = -1            
            vid_axe = -1
            vid_metamer = -1
            vid_node = -1
            vid_elt = -1
            vid_topstem_node = -1
            vid_topstem_element = -1
            metamers = []
            nodes = []
            elts = []
            
        # Add axis
        if axe != prev_axe:
            label = 'axe' + str(axe)           
            if axe == 0:
                vid_axe = g.add_component(vid_plant,edge_type='/',label=label)
            else:
                vid_axe = g.add_child(vid_axe,edge_type='+',label=label)

        # Add metamer
        assert num_metamer > 0
        # args are added to metamers only if metamer_factory is none, otherwise compute metamer components
        args = properties_from_dict(dp,i,exclude=topology)
        components = []
        if metamer_factory:
            if leaf_db:
                xysr = leaf_db[str(args['LcType'])][args['LcIndex']]
            else:
                xysr = {'LcType':args.get('LcType'), 'LcIndex': args.get('LcIndex')}
            components = metamer_factory(Lsect = leaf_sectors, xysr_shape = xysr, **args)
            args={}
        #
        label = 'metamer'+str(num_metamer)
        new_metamer = g.add_component(vid_axe, edge_type='/', label = label, **args)
        if axe==0 and num_metamer==1:
            vid_metamer = new_metamer
        elif num_metamer == 1:
            # add the edge with the bearing metamer on main stem
            vid_metamer = metamers[axe-1]
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
            if axe==0 and num_metamer==1: #root of main stem
                vid_node = new_node
                vid_elt =  new_elt                   
            elif num_metamer == 1: # root of tiller                   
                vid_node = nodes[axe - 1]
                vid_node = g.add_child(vid_node, child = new_node, edge_type='+')
                vid_elt = elts[axe - 1]
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
        if axe == 0 :
            metamers.append(vid_metamer)
            if len(components) > 0:
                nodes.append(vid_topstem_node)
                elts.append(vid_topstem_element)
        prev_plant = plant
        prev_axe = axe
    
    return fat_mtg(g)

