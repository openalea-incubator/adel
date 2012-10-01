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

def internode_elements(Egreen, Esen, diam):
    """ returns parameters of internode elements.
    Fisrt elements is for green visible part of the internode
    Second elements is for senesced visible part of the internode
    """
    green = {'label': 'StemElement', 'length': Egreen, 'tissue_type': 'green', 'diameter' :diam}
    sen = {'label': 'StemElement', 'length': Esen, 'tissue_type': 'senescent', 'diameter' :diam}
    return [green,sen]
    
def sheath_elements(Ggreen, Gsen, diam):
    """ returns parameters of sheath elements.
    Fisrt elements is for green visible part of the sheath
    Second elements is for senesced visible part of the sheath
    """
    green = {'label': 'StemElement', 'length': Ggreen, 'tissue_type': 'green', 'diameter' :diam}
    sen = {'label': 'StemElement', 'length': Gsen, 'tissue_type': 'senescent', 'diameter' :diam}
    return [green,sen]
    
def blade_elements(sectors, Lgreen, Lsen, Laz, Ll, Lw, LcType, LcIndex):
    """ return parameters of blade elements """
    #rolled = {'label': 'StemElement', 'length': Ggreen, 'tissue_type': 'green', 'diameter' :diam}
    ds = 1. / sectors
    Lv = None
    srlim = 1
    if Lv and Lgreen and Lsen:
        Lv = Lgreen + Lsen
        srlim = 1.- Lsen / Lv
    st = ds
    Laz = Laz
    elements=[]
    for isect in range(sectors):
        # create one green and/or one senescent sector
        srb_green = st - ds
        srt_green = min(st, max(srb_green,srlim))
        srb_sen = srt_green
        srt_sen= st
        green = {'label': 'LeafElement', 'length': Lv, 
                'Ll': Ll, 'Lw': Lw, 'LcType': LcType, 'LcIndex': LcIndex,
                'srb': srb_green, 'srt': srt_green, 'green': True}
        sen = {'label': 'LeafElement', 'length': Lv, 
                'Ll': Ll, 'Lw': Lw, 'LcType': LcType, 'LcIndex': LcIndex,
                'srb': srb_sen, 'srt': srt_sen, 'green': False} 
        elements.extend([green,sen])
        st += ds
    return elements
    
def wheat_metamer(Lv=None, Lsen=None, Ll=None, Lw=None, LcType=None, LcIndex=None, Linc=None, Laz=None, Gl=None, Gv=None, Gsen=None, Gd=None, Ginc=None, El=None, Ev=None, Esen=None, Ed=None, Einc=None, *args):
    """ Contructs metamer elements (nodes ?) for adel wheat from parameters returned by setAdel/RunAdel.
    Parameters are : 
      - Lv : visible (emerged) length of blades (green + senesced)
       - Lsen : length of the (emerged) senescent part of the blade
       - Ll : Mature Length of the blade used to compute blade shape (of which the lv first part is exposed)
       - Lw : Maximal width of the mature blade used to compute blade shape
       - LcType : index of a blade shape class
       - LcIndex : index of a blade shape within the class 
       - Linc : relative? inclination of the leaf blade at the top of leaf sheath (deg ?)  
       - Laz : relative ? azimuth of the leaf
        - Gl : Mature length of the sheath
        - Gv : emerged length of the sheath
        - Gsen : emerged senescent length of the sheath
        - Gd : apparent diameter of the sheath
        - Ginc : relative inclination of the sheath
        - El : Mature length of the internode
        - Ev : emerged lenbgthof the internode
        - Esen : emerge senescent length of the internode
        - Ed: diameter of the internode
        - Einc : relative inclination of the internode

    """
    
    Egreen, Ggreen, Lgreen = None, None, None
    
    if Ev and Esen:
        Egreen = Ev - Esen if Ev - Esen > 0. else 0.
        Esen = Ev - Egreen
    if Gv and Gsen:
        Ggreen = Gv - Gsen if Gv - Gsen > 0. else 0.
        Gsen = Gv - Ggreen
    if Lv and Lsen:
        Lgreen = Lv - Lsen if Lv - Lsen > 0. else 0.
        Lsen = Lv - Lgreen
    
   
    modules = [
        {'label': 'internode',
        'green_visible_length': Egreen,
        'senesced_visible_length': Esen,
        'mature_length': El,
        'diameter' : Ed, 
        'inclination' : Einc,
        'elements' : internode_elements(Egreen, Esen, Ed)}, 
        {'label': 'sheath',
        'green_visible_length': Ggreen,
        'senesced_visible_length': Gsen,
        'mature_length': Gl,
        'diameter' : Gd,
        'inclination' : Ginc,
        'azimuth' : Laz,
        'elements': sheath_elements(Ggreen,Gsen,Gd)}, 
        {'label': 'blade', 
        'green_visible_length': Lgreen,
        'senesced_visible_length': Lsen,
        'mature_length': Ll,
        'max_width' : Lw,
        'inclination' : Linc,
        'LcType': LcType,
        'LcIndex' : LcIndex,
        'elements': blade_elements(1,Lgreen, Lsen, Laz, Ll, Lw, LcType, LcIndex)} 
    ]
    
    return modules
    
    
    
def mtg_factory(parameters, metamer_factory=None, leaf_sectors=1, topology = ['plant','axe','numphy']):
    """ Construct a MTG from a dictionary of parameters.

    The dictionary contains the parameters of all metamers in the stand (topology + properties).
    metamer_factory is a function that build metamer properties and metamer elements from parameters. It is called with all metamer properties, 
    and should return a dict of elements. An element is a dict of properties to be attached to that element.
    Valid elements name  are : 'metamer', 'internode_xxx', 'sheath_xxx', and 'blade_xxx', with xxx being the external appearance ('hidden', 'green' or 'senescent')
    Elements missing in metamer_factory output (or if metammer factory is none), are created without properties
    Sector is an integer giving the number of LeafElements per Leaf blade
    add_hidden adds hidden parts of elements
    topology is the list of key names used in parameters for plant number, axe numer and metamer number

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
            components = metamer_factory(**args)
            args={}
        #
        label = 'metamer'+str(num_metamer)
        if axe==0 and num_metamer==1:
            vid_metamer = g.add_component(vid_axe, edge_type='/', label = label, **args)
        elif num_metamer == 1:
            # create as the first component of the axilary axis
            new_metamer = g.add_component(vid_axe, label=label, edge_type='/', **args)
            # add the edge with the bearing metamer on main stem
            vid_metamer = metamers[axe-1]
            vid_metamer =  g.add_child(vid_metamer, child=new_metamer, edge_type='+')
            #vid_node = nodes[axe-1]
        else:
            vid_metamer = g.add_child(vid_metamer, label=label, edge_type='<', **args)

        # add metamer components, if any           
        if len(components) > 0:
            def get_component(components, index):
                component = components[index]
                elements = component['elements']
                properties = dict(component)
                del properties['elements']
                return elements, properties
            # deals with first component and first element
            elements, properties = get_component(components,0)
            element = elements[0]
            new_node = g.add_component(vid_metamer, edge_type='/', **properties)
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
                elements, properties = get_component(components,i)
                edge_type = '<'
                if properties['label'] == 'sheath':
                    edge_type = '+'
                vid_node = g.add_child(vid_node, edge_type=edge_type, **properties)      
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

