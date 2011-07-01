# -*- coding: utf-8 -*-

"""
Methods on the mtg for adel programs.
Author: Christophe Pradal
License: CeCILL
Date: 08/07/2008
"""

import csv

from openalea.mtg import MTG, fat_mtg

try:
    from openalea.mtg.traversal import *
except:
    pass

from openalea.mtg.io import read_lsystem_string
from openalea.mtg.algo import union

import numpy
try:
    from openalea.plantgl.all import (Scene,Translation, Vector3, Geometry,
                                AxisRotation, AxisRotated,  Transform4, BaseOrientation, 
                                Shape, Material, Color3, PglTurtle, Mesh, Translated) 
except:
    Material = tuple
    Color3 = tuple
    pass
import random
from math import pi

def to_aggregation_table(g):
    """
    Convert a mtg  `g` to an aggregation table.
    Returns a multi-line string with one line per triangle.
    """

     
    # 1. classe at each scale
    label = g.property('label')
    index = g.property('index')
    geometry = g.property('geometry')
    can_label = g.property("can_label")

    symbols = dict(zip(label.itervalues(), label.iterkeys()))
    for k, v in symbols.iteritems():
        symbols[k] = g.scale(v)

    l = symbols.items()
    l.sort(cmp= lambda x, y: cmp(x[1], y[1]))

    # scales
    header = "Plant Axe Metamer StemElement LeafElement Type"
    header = header.split()

    # 2. iterate on the geometry at the last scale
    assert g.max_scale() == 4

    root_axe = list(g.roots(scale=2))[0]
    root_metamer = list(g.roots(scale=3))[0]
    root_elt = list(g.roots(scale=4))[0]

    # compute the number of triangles

    nb_lines = sum( mesh.indexListSize() for mesh in geometry.itervalues() if mesh)
    lines = numpy.zeros((nb_lines, 6), dtype=int)

    # compute relative index for metamer in axe and element in metamer
    local_index = {}
    #determine the number of element
    nb_stem_elements = {}
    for root_axe in g.roots(scale=2):
        for axe_id in pre_order(g,root_axe):
            for i, mid in enumerate(g.components(axe_id)):
                local_index[mid] = i+1

    for root_metamer in g.roots(scale=3):
        for mid in pre_order(g,root_metamer):
            stem_index = 0
            leaf_index = 0
            for eid in g.components(mid):
                if 'stem' in label[eid].lower():
                    stem_index += 1
                    local_index[eid] = stem_index
                else:
                    leaf_index += 1
                    local_index[eid] = leaf_index
            
            nb_stem_elements[mid] = stem_index

    i = 0
    for root_elt in g.roots(scale=4):
        for vid in pre_order(g, root_elt):
            metamer_id = g.complex(vid)
            axe_id = g.complex(metamer_id)
            plant_id = g.complex(axe_id)

            plant_index = index[plant_id]
            axe_index = index[axe_id]
            metamer_index = local_index[metamer_id]
            element_index = local_index[vid]

            stem_index = 0
            leaf_index = 0

            # element
            is_stem = 'stem' in label[vid].lower()
            if is_stem:
                stem_index = element_index
            else:
                leaf_index = element_index
            #determine tissue type
            lab = can_label[vid]
            if lab.is_soil():
                ttype = 0
            elif lab.is_leaf():
                ttype = 1#lamina
            elif lab.is_stem():
                if nb_stem_elements[metamer_id] == 1:
                    ttype = 2 if lab.optical_id <= 2 else 5#2 = sheath, 5 = ear
                else: 
                    if element_index == 1:
                        ttype = 3 if lab.optical_id <= 2 else 4#internode(4) or peduncle(4)
                    else:
                        ttype = 2 if lab.optical_id <= 2 else 5#sheath or ear
            else:
                ttype = -1 #unknown
             

            
            indices = (plant_index, axe_index, metamer_index, stem_index, leaf_index, ttype)
            geom = geometry.get(vid)
            if geom:
                n = geom.indexListSize()
                for j in range(n):
                    lines[i] = indices
                    i += 1

    lines = lines.transpose()
    cols =[c.transpose() for c in lines]
    table = dict(zip(header, cols))
    return table
        
def to_plantgl(g, 
               leaf_material = None,
               stem_material = None,
               soil_material = None):
    """
    Returns a plantgl scene from an mtg.
    """
    if leaf_material is None:
        leaf_material = Material(Color3(0,180,0))
    if stem_material is None:
        stem_material = Material(Color3(0,130,0))
    if soil_material is None:
        soil_material = Material(Color3(170, 85,0))

    geometries = g.property('geometry')
    labels = g.property('can_label')
    scene = Scene()


    def geom2shape(vid, mesh, scene):
        shape = None
        if isinstance(mesh, list):
            for m in mesh:
                geom2shape(vid, m, scene)
            return
        if mesh is None:
            return
        if isinstance(mesh, Shape):
            shape = mesh
            mesh = mesh.geometry
        label = labels.get(vid)
        if not label:
            if not shape:
                shape = Shape(mesh)
        elif label.is_soil():
            shape = Shape(mesh, soil_material)
        elif label.is_stem() and label.optical_id <= 1:
            shape = Shape(mesh, stem_material)
        elif label.is_stem() and label.optical_id > 1:
            shape = Shape(mesh, soil_material)
        elif label.is_leaf() and label.optical_id <= 1:
            shape = Shape(mesh, leaf_material)
        elif label.is_leaf() and label.optical_id > 1:
            shape = Shape(mesh, soil_material)
        shape.id = vid
        scene.add(shape)

    for vid, mesh in geometries.iteritems():
        geom2shape(vid, mesh, scene)
    return scene


def _line(ind, pts, label):
    s = "p 1 %s 3 %s"%(str(label), ' '.join('%.6f'%x for i in ind for x in pts[i]))
    return s

def to_canestra(g):
    """
    Return a string representing a canestra file.
    """
    geometry = g.property('geometry')
    can_label = g.property("can_label")

    begin = '# File generated by OpenAlea.Adel program'

    lines = [begin]
    max_scale = g.max_scale()

    for root_elt in g.roots(scale=max_scale):
        for vid in pre_order(g, root_elt):
            mesh = geometry.get(vid)
            if not mesh:
                continue
            pts = numpy.array(mesh.pointList, ndmin=2)
            indices = numpy.array(mesh.indexList, ndmin=2)
            label = can_label[vid]
            lines.extend([_line(ind, pts, label) for ind in indices])
    lines.append('')
    return '\n'.join(lines)

def planter(g, distribution, random_seed=0):
    """
    Transform a set of plants with given point distributions.

    :Parameters: 
        - g: MTG
        - distribution: a list of 2D positions
        - random_seed: add a random rotation to each plant if the value is positive.
          In this case, random_seed is used as a seed.
    """
    #assert g.nb_vertices(scale=1) == len(distribution)
    geometry = g.property('geometry')
    
    # store the previous plant translation to not copy each time all the geometry.
    if '_plant_translation' not in g.properties():
        g.add_property('_plant_translation')
    translations = g.property('_plant_translation')
    

    if random_seed > 0:
        random.seed(random_seed)

    def pt2transfo(pt, previous_pt):
        matrix = Translation(pt).getMatrix()
        rotation = 0
        if random_seed > 0:
            rotation = random.random()*pi
            r = AxisRotation((0,0,1), rotation).getMatrix()
            matrix = matrix * r * Translation(previous_pt).getMatrix()
        return Transform4(matrix), pt, rotation

    def transform_geom(geom, transfo, translation, rotation):
        if isinstance(geom, list):
            geom = [transform_geom(g,transfo, translation, rotation) for g in geom]
        elif isinstance(geom, Mesh):
            geom = geom.transform(transfo) if geom else geom
        elif isinstance(geom, Geometry):
            geom = Translated(translation, AxisRotated((0,0,1),rotation, geom))
        elif isinstance(geom, Shape):
            geom = Shape(Translated(translation, AxisRotated((0,0,1),rotation, geom.geometry)))
        return geom

    plants = g.vertices(scale=1)
    plants = list(plants)[:len(distribution)]

    max_scale = g.max_scale()

    for i, root_elt in enumerate(plants):
        previous_translation = translations.get(root_elt,(0,0,0))
        translations[root_elt] = distribution[i]

        #displacement = Vector3(distribution[i])-Vector3(previous_translation)

        transfo, translation, rotation = pt2transfo(Vector3(distribution[i]), -Vector3(previous_translation))
        l = (vid for vid in g.vertices(scale=max_scale) if g.complex_at_scale(vid, scale=1) == root_elt) 
        #for vid in g.components_at_scale(root_elt, 4):
        for vid in l:
            geom = geometry.get(vid)
            geom = transform_geom(geom, transfo, translation, rotation)

            if geom:
                geometry[vid] = geom
                can_label = g.property('can_label').get(vid)
                if can_label:
                    can_label.plant_id = i

def duplicate(g, n=1):
    if g.nb_vertices(scale=1) == n:
        return g
    g1 = g.sub_mtg(g.root)
    for i in range(n-1):
        g1 = union(g1,g)
    return g1

class CanMTG(MTG):
    def __init__(self, functions, axial_string):
        MTG.__init__(self)
        symbols = {'newPlant' : 1, 'newAxe' : 2, 'newMetamer' :3, 'StemElement':4, 'LeafElement':4}
        g = read_lsystem_string(axial_string, symbols, functions, self)
        self.symbols = symbols

CanMTG.planter = planter
CanMTG.to_plantgl = to_plantgl
CanMTG.to_canestra = to_canestra
CanMTG.to_aggregation_table = to_aggregation_table

MTG.planter = planter
MTG.to_plantgl = to_plantgl
MTG.to_canestra = to_canestra
MTG.to_aggregation_table = to_aggregation_table

def convert(v, undef='NA'):
    res = v
    try:
        res = int(res)
    except ValueError:
        try:
            res = float(res)
        except ValueError:
            if res == undef:
                res = None
    return res

def properties(d, exclude = []):
    res = {}
    for k, v in d.iteritems():
        if k in exclude:
            continue
        v = convert(v)
        if v is not None:
            res[k] = v
    return res

def properties_from_dict(d,index, exclude = []):
    res = {}
    for k in d:
        if k in exclude:
            continue
        res[k] = d[k][index]
    return res

def topological_table_to_mtg(csv_file, epsilon=1e-6):
    f = open(csv_file)
    l=f.readline()

    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(l)
    # Go from start of file
    f.seek(0)

    reader = csv.DictReader(f, dialect=dialect)

    # Build the MTG
    g = MTG()
    topology = ['plant', 'axe', 'numphy']

    # buffers
    prev_plant = 0
    prev_axe = -1
    prev_metamer = -1
    vid_plant = -1
    vid_axe = -1
    vid_metamer = -1
    vid_node = -1

    # store the metamer vids of the axis 0 
    metamers = []
    nodes = []
    edge_type = '<'
    for d in reader:
        plant, axe, num_metamer = [int(convert(d.get(x),undef=None)) for x in topology]
        print 'adding vertices for plant: %d, axe:%d, metamer: %d'%(plant, axe, num_metamer)
        # Add new plant
        if plant != prev_plant:
            label = 'plant'+str(plant)
            vid_plant = g.add_component(g.root, edge_type='/',label=label)
            vid_axe = -1
            vid_metamer = -1
            vid_node = -1
            prev_axe = -1
        #??
        if num_metamer < prev_metamer:
            prev_axe = -1
            prev_metamer = -1
            
        # Add an axis
        if axe != prev_axe:
            label = 'axe'+str(axe)
           
            
            if axe == 0:
                vid_axe = g.add_component(vid_plant,edge_type='/',label=label)
                vid_node = -1
            else:
                #args['edge_type'] = '+'
                edge_type = '+'
                vid_axe = g.add_child(vid_axe, edge_type=edge_type, label=label)
                #vid, vid_axe = g.add_child_and_complex(metamers[axe], complex=vid_axe, **args)

        # Add metamers
        args = properties(d, exclude=topology)

        assert num_metamer > 0
        label = 'metamer'+str(num_metamer)
        new_axe = True

        if axe==0 and num_metamer==1:
            metamers=[]
            edge_type = '/'
            vid_metamer = g.add_component(vid_axe, edge_type='/', label = label, **args)
        elif num_metamer == 1:
            # Add the first element connected to previous metamer on the previous axis
            edge_type = '+'
            vid_metamer, vid_axe =  g.add_child_and_complex(metamers[axe-1], complex=vid_axe, edge_type='+',label=label, **args)
            vid_node = nodes[axe-1]
        else:
            edge_type = '<'
            vid_metamer = g.add_child(vid_metamer, label=label, edge_type='<',**args)
            new_axe = False

        if axe == 0:
            metamers.append(vid_metamer)

        # Add internode, sheath, and lamina
        nb_internodes = 0
        nb_sheath = 0
        nb_leaf = 0

        # Internode
        if args['Ev'] > 0.:
            if args['Esen'] > args['Ev']:
                if vid_node == -1:
                    vid_node= g.add_component(vid_metamer, label='Esen', edge_type='/',length=args['Ev'],po=args['Epos'], diam=args['Ed'] )
                    assert edge_type == '/'
                else:
                    vid_node, vid_metamer= g.add_child_and_complex(vid_node, complex=vid_metamer, label='Esen', edge_type=edge_type,length=args['Ev'],po=args['Epos'], diam=args['Ed'] )
            else:
                if vid_node == -1:
                    vid_node= g.add_component(vid_metamer, label='Egreen', edge_type='/',length=args['Ev']-args['Esen'],po=args['Epo'], diam=args['Ed'] )
                    assert edge_type == '/'
                else:
                    vid_node, vid_metamer= g.add_child_and_complex(vid_node, complex=vid_metamer, label='Egreen', edge_type=edge_type,length=args['Ev']-args['Esen'],po=args['Epo'], diam=args['Ed'] )
                vid_node = g.add_child(vid_node, label='Esen', edge_type='<',length=args['Esen'],po=args['Epos'], diam=args['Ed'])
        else:
            if new_axe and args['Gv'] == 0.:
                if vid_node == -1:
                    vid_node= g.add_component(vid_metamer, label='Egreen', edge_type='/',length=0.,po=args['Epo'], diam=args['Ed'] )
                else:
                    vid_node, vid_metamer= g.add_child_and_complex(vid_node, complex=vid_metamer, label='Egreen', edge_type=edge_type,length=0.,po=args['Epo'], diam=args['Ed'] )
        # Sheath
        if args['Gv'] > 0.:
            if args['Gsen'] > args['Gv']:
                vid_node= g.add_child(vid_node, label='Gsen', edge_type='<',length=args['Gv'],po=args['Gpos'], diam=args['Gd'] )
            else:
                vid_node = g.add_child(vid_node, label='Ggreen', edge_type='<',length=args['Gv']-args['Gsen'],po=args['Gpo'], diam=args['Gd'] )
                vid_node = g.add_child(vid_node, label='Gsen', edge_type='<',length=args['Gsen'],po=args['Gpos'], diam=args['Gd'])

        if axe == 0:
            nodes.append(vid_node)

        # Laminae
        if args['Lv'] > 0.:
            if args['Lsen'] > args['Lv']:
                l_node= g.add_child(vid_node, label='Lsen', edge_type='+',length=args['Lv'],po=args['Lpos'],
                                      Ll=args['Ll'], Lw= args['Lw'], LcType=args['LcType'], LcIndex=args['LcIndex'], Linc=args['Linc'], Laz=args['Laz'],
                                      srb=0, srt=1)
            else:
                l_node = g.add_child(vid_node, label='Lgreen', edge_type='+',length=args['Lv']-args['Lsen'],po=args['Lpo'],
                                       Ll=args['Ll'], Lw= args['Lw'], LcType=args['LcType'], LcIndex=args['LcIndex'], Linc=args['Linc'], Laz=args['Laz'],
                                      srb=0, srt=1-(args['Lsen']/args['Lv']))
                l_node = g.add_child(l_node, label='Lsen', edge_type='<',length=args['Lsen'],po=args['Lpos'],
                                       Ll=args['Ll'], Lw= args['Lw'], LcType=args['LcType'], LcIndex=args['LcIndex'], Linc=args['Linc'], Laz=args['Laz'],
                                      srt=1, srb=1-(args['Lsen']/args['Lv']))

        
        prev_plant = plant
        prev_axe = axe
        prev_metamer = num_metamer

        


    f.close()

    return fat_mtg(g)

def internode(g, vid_axe, vid_metamer, prev_node, props, epsilon):
    if props['Ev'] < epsilon:
        return prev_node

def mtg_factory(params):
    """ Construct a MTG from a dictionary.

    The dictionary contains the parameters and a list of elements.
    The keys of params are:
        - plant: indx of plant
        - axe: 
        - numphy
        - Lv
        - Ll
        - Lsen
        - Lw
        - LcType
        - LcIndex
        - Linc
        - Laz
        - Lpo
        - Lpos
        - Gl
        - Gv
        - Gsen
        - Gpo
        - Gpos
        - Gd
        - Ginc
        - El
        - Ev
        - Esen
        - Ed
        - Einc
        - Epo
        - Epos

    :TODO: 
        * add length and final_length (DONE)
        * fix naming convention for Linc: relative inclination (DONE)
        * Add a scale to define the tissue types (ear, sheath, laminae, gain)
        * diam and final_diam (resp. width)
        * function reset length
        * function phenology(g, table) -> dynamic parameters (start_thermaltime, end_thermaltime)
        * function growth_thermaltime(g, tt, d_tt): tt=thermaltime du couvert
        * function growth_thermaltime(g, tt, d_tt, stress factor)
        * stress_factor: offre/demande
            - demand = :math:`D=\int_{tt}^{tt+dtt}{S(x)dx}*\rho_s+\int_{tt}^{tt+dtt}{V(x)dx}*\rho_v`
            - offre : :math:`sum{E_abs}*\eps_b`
            => ds = ds_predit* stress_factor
        * give the area to the leaf model
        * update properties
    """

    if not _check_adel_parameters(params):
        raise ValueError('Adel parameters are invalid')

    g = MTG()
    topology = ['plant', 'axe', 'numphy']

    # buffers
    prev_plant = 0
    prev_axe = -1
    prev_metamer = -1
    vid_plant = -1
    vid_axe = -1
    vid_metamer = -1
    vid_node = -1

    # store the metamer vids of the axis 0 
    metamers = []
    nodes = []
    edge_type = '<'

    dp = params
    nrow = len(params['plant'])
    for i in range(nrow):

        #plant, axe, num_metamer = [int(convert(d.get(x),undef=None)) for x in topology]
        plant = int(dp['plant'][i])
        axe = int(dp['axe'][i])
        num_metamer = int(dp['numphy'][i])

        #plant, axe, num_metamer = [int(convert(d.get(x),undef=None)) for x in topology]
        #print 'plant: %d, axe:%d, nb_metamers: %d'%(plant, axe, num_metamer)
        # Add new plant
        if plant != prev_plant:
            label = 'plant'+str(plant)
            vid_plant = g.add_component(g.root, label=label, edge_type='/')
            vid_axe = -1
            vid_metamer = -1
            vid_node = -1
            prev_axe = -1

        if num_metamer < prev_metamer:
            prev_axe = -1
            prev_metamer = -1
            
        # Add an axis
        if axe != prev_axe:
            label = 'axe'+str(axe)
           
            
            if axe == 0:
                vid_axe = g.add_component(vid_plant,edge_type='/',label=label)
                vid_node = -1
            else:
                #args['edge_type'] = '+'
                edge_type = '+'
                assert g.parent(vid_plant) is None
                new_axe = g.add_child(vid_axe,edge_type=edge_type,label=label)
                print vid_axe, new_axe, vid_plant
                assert g.parent(vid_plant) is None
                #vid, vid_axe = g.add_child_and_complex(metamers[axe], complex=vid_axe, **args)

        # Add metamers
        args = properties_from_dict(dp,i,exclude=topology)

        assert num_metamer > 0
        label = 'metamer'+str(num_metamer)

        if axe==0 and num_metamer==1:
            metamers=[]
            edge_type = '/'
            vid_metamer = g.add_component(vid_axe, edge_type='/', label = label, **args)
        elif num_metamer == 1:
            # Add the first element connected to previous metamer on the previous axis
            edge_type = '+'
            vid_metamer = g.add_component(vid_axe, label=label, **args)
            vid_metamer =  g.add_child(metamers[axe-1], child=vid_metamer, edge_type='+')
            vid_node = nodes[axe-1]
        else:
            edge_type = '<'
            vid_metamer = g.add_child(vid_metamer, label=label, edge_type='<',**args)

        if axe == 0:
            metamers.append(vid_metamer)

        # Add internode, sheath, and lamina

        # Visible Internode
        if args['Ev'] > 0.:
            tissue = 'internode'

            # Several cases:
            if args['Esen'] > args['Ev']:
                #  - Only one senescent element
                if vid_node == -1:
                    # first metamer on the axis
                    vid_node= g.add_component(vid_metamer, label='StemElement', edge_type='/',
                        length=args['Ev'],final_length=args['El'], po=args['Epos'], diam=args['Ed'], tissue=tissue, sen=True, inclination=args['Einc'] )
                    assert edge_type == '/'
                else:
                    # first element on the metamer which has a parent metamer on the axis
                    new_node = g.add_component(vid_metamer)
                    vid_node= g.add_child(vid_node,child=new_node, 
                        label='StemElement', edge_type=edge_type,
                        length=args['Ev'], final_length=args['El'], po=args['Epos'], diam=args['Ed'], tissue=tissue, sen=True,inclination=args['Einc'] )
            else:
                #  - Green element and perhaps a senescent one
                if vid_node == -1:
                    # first metamer on the axis
                    vid_node= g.add_component(vid_metamer, label='StemElement', edge_type='/',
                        length=args['Ev']-args['Esen'],final_length=args['El'],po=args['Epo'], diam=args['Ed'], tissue=tissue, inclination=args['Einc'] )
                    assert edge_type == '/'
                else:
                    # first element on the metamer which has a parent metamer on the axis
                    new_node = g.add_component(vid_metamer)
                    vid_node= g.add_child(vid_node,child=new_node, 
                            label='StemElement', edge_type=edge_type,
                            length=args['Ev']-args['Esen'],final_length=args['El'],po=args['Epo'], diam=args['Ed'], tissue=tissue, inclination=args['Einc'] )
                # Add a senescent element
                if args['Esen'] > 0.:
                    vid_node = g.add_child(vid_node, label='StemElement', edge_type='<',
                    length=args['Esen'],final_length=args['El'],po=args['Epos'], diam=args['Ed'], tissue=tissue, sen=True, inclination=0.)
        else:
            tissue = 'internode'
            if args['Gv'] <= 0.:
                if vid_node == -1:
                    # first metamer on the axis
                    vid_node= g.add_component(vid_metamer, label='StemElement', edge_type='/',
                        length=args['Ev'],final_length=args['El'], po=args['Epos'], diam=args['Ed'], tissue=tissue, sen=True, inclination=args['Einc'] )
                    assert edge_type == '/'
                else:
                    # first element on the metamer which has a parent metamer on the axis
                    new_node = g.add_component(vid_metamer)
                    vid_node= g.add_child(vid_node,child=new_node, 
                        label='StemElement', edge_type=edge_type,
                        length=args['Ev'], final_length=args['El'], po=args['Epos'], diam=args['Ed'], tissue=tissue, sen=True,inclination=args['Einc'] )
               
        # Sheath
        # Same logic that described previously.
        if args['Gv'] > 0.:
            tissue = 'sheath'
            if args['Gsen'] > args['Gv']:
                if vid_node == -1:
                    vid_node= g.add_component(vid_metamer, label='StemElement', edge_type='/',
                        length=args['Gv'],final_length=args['Gl'],po=args['Gpos'], diam=args['Gd'], tissue=tissue, sen=True, inclination=args['Ginc'])
                else:
                    new_node = g.add_component(vid_metamer)
                    vid_node= g.add_child(vid_node, child=new_node, label='StemElement', edge_type='<',
                    length=args['Gv'],final_length=args['Gl'],po=args['Gpos'], diam=args['Gd'], tissue=tissue, sen=True, inclination=args['Ginc'] )
            else:
                if vid_node == -1:
                    vid_node = g.add_component(vid_metamer, label='StemElement', edge_type='/',
                    length=args['Gv']-args['Gsen'],final_length=args['Gl'],po=args['Gpo'], diam=args['Gd'], tissue=tissue, inclination=args['Ginc'] )
                else:
                    new_node = g.add_component(vid_metamer)
                    vid_node= g.add_child(vid_node, child=new_node, label='StemElement', edge_type='<',
                    length=args['Gv']-args['Gsen'],final_length=args['Gl'],po=args['Gpo'], diam=args['Gd'], tissue=tissue, inclination=args['Ginc'] )
                if args['Gsen'] > 0.:
                    vid_node= g.add_child(vid_node, label='StemElement', edge_type='<',
                    length=args['Gsen'],final_length=args['Gl'],po=args['Gpos'], diam=args['Gd'], sen=True, tissue=tissue, inclination=0.)

        if axe == 0:
            nodes.append(vid_node)

        # Lamina
        if args['Lv'] > 0.:
            tissue='lamina'
            # Only one senescent LeafElement
            if args['Lsen'] > args['Lv']:
                l_node= g.add_child(vid_node, label='LeafElement', edge_type='+',length=args['Lv'],final_length=args['Ll'],
                    po=args['Lpos'], Ll=args['Ll'], Lw=args['Lw'], 
                    LcType=args['LcType'], LcIndex=args['LcIndex'], inclination=0.,Linc=args['Linc'], 
                    Laz=args['Laz'], srb=0, srt=1, tissue=tissue, sen=True)
            else:
                # one green LeafElement followed eventually by one senescent
                l_node = g.add_child(vid_node, label='LeafElement', edge_type='+',
                    length=args['Lv']-args['Lsen'],final_length=args['Ll'],po=args['Lpo'],
                    Ll=args['Ll'], Lw= args['Lw'], LcType=args['LcType'], 
                    LcIndex=args['LcIndex'], inclination=0., Laz=args['Laz'],Linc=args['Linc'], 
                    srb=0, srt=1-(args['Lsen']/args['Lv']), tissue=tissue)
                if args['Lsen'] > 0.:
                    l_node = g.add_child(l_node, label='LeafElement', edge_type='<',
                        length=args['Lsen'],final_length=args['Ll'],po=args['Lpos'],
                        Ll=args['Ll'], Lw= args['Lw'], LcType=args['LcType'], 
                        LcIndex=args['LcIndex'], inclination=0., Laz=args['Laz'],Linc=args['Linc'], 
                        srt=1, srb=1-(args['Lsen']/args['Lv']), sen=True, tissue=tissue)

        
        prev_plant = plant
        prev_axe = axe
        prev_metamer = num_metamer

    
    return fat_mtg(g)


def _check_adel_parameters( params):
    return True

def compute_element(n, symbols):
    leaf = symbols.get('LeafElement')
    stem = symbols.get('StemElement')

    leaf_rank = int(n.complex().index())
    optical_species = int(n.po)
    final_length = n.final_length
    length = n.length
    s_base = n.srb
    s_top = n.srt
    seed = n.LcIndex
    #relative leaf inclination
    linc = n.Linc
    
    element = {} 
    if n.label.startswith('L'):
        radius_max = n.Lw
        element = leaf(optical_species, 
                    final_length, 
                    length, 
                    radius_max, 
                    s_base, 
                    s_top, 
                    leaf_rank,
                    seed,
                    linc) 
    else:
        diameter_base = n.parent().diam if (n.parent() and n.parent().diam > 0.) else n.diam
        diameter_top = n.diam
        element = stem( optical_species, length, diameter_base, diameter_top)

    can_label =  element['label']
    if can_label:
        can_label.elt_id = leaf_rank
        plant_node = n.complex_at_scale(scale=1)
        can_label.plant_id = plant_node.index()

    return element['geometry'], can_label

def transform(turtle, mesh):
        x = turtle.getUp()
        z = turtle.getHeading()

        bo = BaseOrientation(x, z^x)
        matrix = Transform4(bo.getMatrix())
        matrix.translate(turtle.getPosition())
        mesh = mesh.transform(matrix)
        return mesh





def mtg_turtle(g, symbols):
    ''' Compute the geometry on each node of the MTG using Turtle geometry. '''

    from openalea.mtg import turtle

    plants = g.component_roots_at_scale(g.root, scale=1)
    nplants = g.nb_vertices(scale=1)

    gt = MTG()

    def adel_visitor(g, v, turtle):
        # 1. retriev the node
        n = g.node(v)
        angle = float(n.Laz) if n.Laz else 0.
        turtle.rollL(angle)
        if g.edge_type(v) == '+':
            if not n.label.startswith('L'):
                angle = n.Ginc or n.Einc
                angle = float(angle) if angle is not None else 0.
                turtle.up(angle)

        # 2. Compute the geometric symbol
        mesh, can_label = compute_element(n, symbols)
        if mesh:
            n.geometry = transform(turtle, mesh)
            n.can_label = can_label

        # 3. Update the turtle
        turtle.setId(v)
        if n.label.startswith('S'):
            turtle.f(n.length)
        # Get the azimuth angle

    for plant in plants:
       gplant = g.sub_mtg(plant)
       scene = turtle.TurtleFrame(gplant,visitor=adel_visitor)
       gt = union(gplant,gt)
       
    return gt
    

def mtg_turtle_time(g, symbols, time):
    ''' Compute the geometry on each node of the MTG using Turtle geometry. '''

    from openalea.mtg import turtle

    g.properties()['geometry'] = {}
    g.properties()['_plant_translation'] = {}

    def compute_element(n, symbols, time):
        leaf = symbols.get('LeafElement')
        stem = symbols.get('StemElement')

        leaf_rank = int(n.complex().index())
        optical_species = int(n.po)
        final_length = n.final_length
        try :
            length = final_length * (time - n.start_tt) / (n.end_tt - n.start_tt) if n.end_tt and time < n.end_tt else n.length
        except:
            length = n.length
        s_base = n.srb
        s_top = n.srt
        seed = n.LcIndex
        #relative leaf inclination
        linc = n.Linc
 
        element = {} 
        if n.label.startswith('L'):
            radius_max = n.Lw
            element = leaf(optical_species, 
                        final_length, 
                        length, 
                        radius_max, 
                        s_base, 
                        s_top, 
                        leaf_rank, seed, linc) 
        else:
            diameter_base = n.parent().diam if (n.parent() and n.parent().diam > 0.) else n.diam
            diameter_top = n.diam
            element = stem( optical_species, length, diameter_base, diameter_top)

        can_label =  element['label']
        if can_label:
            can_label.elt_id = leaf_rank
            plant_node = n.complex_at_scale(scale=1)
            can_label.plant_id = plant_node.index()

        return element['geometry'], can_label

    def adel_visitor(g, v, turtle, time):
        # 1. retriev the node

        n = g.node(v)
        angle = float(n.Laz) if n.Laz else 0.
        turtle.rollL(angle)
        if g.edge_type(v) == '+':
            angle = n.inclination
            angle = float(angle) if angle is not None else 0.
            turtle.up(angle)

        # 2. Compute the geometric symbol
        mesh, can_label = compute_element(n, symbols, time)
        if mesh:
            n.geometry = transform(turtle, mesh)
            n.can_label = can_label

        # 3. Update the turtle
        turtle.setId(v)

        try:
            length = n.length * (time - n.start_tt) / (n.end_tt - n.start_tt) if time < n.end_tt else n.length
        except:
            length = n.length
	if ('Leaf' not in n.label) and (length > 0.):
        	turtle.F(length)
        # Get the azimuth angle
        

    def traverse_with_turtle_time(g, vid, time, visitor=adel_visitor):
        turtle = PglTurtle()
        times = g.property('time')
        def push_turtle(v):
            n = g.node(v)
            try:
                if n.start_tt > time:
                    return False
            except: 
                pass
            if g.edge_type(v) == '+':
                turtle.push()
            return True

        def pop_turtle(v):
            n = g.node(v)
            try:
                if n.start_tt > time:
                    return 
            except: 
                pass
            if g.edge_type(v) == '+':
                turtle.pop()

        visitor(g,vid,turtle,time)
        turtle.push()
	plant_id = g.complex_at_scale(vid, scale=1)
        for v in pre_order2_with_filter(g, vid, None, push_turtle, pop_turtle):
            if v == vid: continue
            visitor(g,v,turtle,time)

        scene = turtle.getScene()
        return g

    for plant_id in g.component_roots_at_scale(g.root, scale=4):
        g = traverse_with_turtle_time(g, plant_id, time)
    return g

def thermal_time(g):
    """ Dummy function to test adel with a thermal time parameter.
    """

    plants = g.vertices(scale=1)
    metamer_scale = g.max_scale()-1
    dtt = 10.

    for plant in plants:
        tt = 0
        v = g.component_roots_at_scale(g.root, scale=metamer_scale).next()
        for metamer in pre_order2(g, v):
            nm = g.node(metamer)
            for node in nm.components():
                node.start_tt = tt
                node.end_tt = tt+dtt
            tt += dtt
    return g

