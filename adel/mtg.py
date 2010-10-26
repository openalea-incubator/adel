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
import numpy
try:
    from openalea.plantgl.all import (Scene,Translation, Vector3, 
                                AxisRotation, Transform4, 
                                Shape, Material, Color3) 
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

    symbols = dict(zip(label.itervalues(), label.iterkeys()))
    for k, v in symbols.iteritems():
        symbols[k] = g.scale(v)

    l = symbols.items()
    l.sort(cmp= lambda x, y: cmp(x[1], y[1]))

    # scales
    header = "Plant Axe Metamer StemElement LeafElement"
    header = header.split()

    # 2. iterate on the geometry at the last scale
    assert g.max_scale() == 4

    root_axe = list(g.roots(scale=2))[0]
    root_metamer = list(g.roots(scale=3))[0]
    root_elt = list(g.roots(scale=4))[0]

    # compute the number of triangles

    nb_lines = sum( mesh.indexListSize() for mesh in geometry.itervalues() if mesh)
    lines = numpy.zeros((nb_lines, 5), dtype=int)

    # compute relative index for metamer in axe and element in metamer
    local_index = {}
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

            indices = (plant_index, axe_index, metamer_index, stem_index, leaf_index)
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

    for vid, mesh in geometries.iteritems():
        if mesh is None:
            continue
        label = labels.get(vid)
        if not label:
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
    
    for root_elt in g.roots(scale=4):
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

def planter(g, distribution):
    """
    Transform a set of plants with given point distributions.
    """
    #assert g.nb_vertices(scale=1) == len(distribution)
    geometry = g.property('geometry')
    
    # store the previous plant translation to not copy each time all the geometry.
    if '_plant_translation' not in g.properties():
        g.add_property('_plant_translation')
    translations = g.property('_plant_translation')
    
    def pt2transfo(pt):
        #r = AxisRotation((0,0,1), random.random()*pi).getMatrix()
        #return Transform4(Translation(pt).getMatrix()*r)
        return Transform4(Translation(pt).getMatrix())

    plants = g.vertices(scale=1)
    plants = list(plants)[:len(distribution)]

    for i, root_elt in enumerate(plants):
        previous_translation = translations.get(root_elt,(0,0,0))
        translations[root_elt] = distribution[i]
        displacement = Vector3(distribution[i])-Vector3(previous_translation)

        transfo = pt2transfo(displacement)
        l = (vid for vid in g.vertices(scale=4) if g.complex_at_scale(vid, scale=1) == root_elt) 
        #for vid in g.components_at_scale(root_elt, 4):
        for vid in l:
            geom = geometry.get(vid)
            if geom:
                geometry[vid] = geom.transform(transfo)
            else:
                geometry[vid] = geom

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
            if vid_node == -1:
                vid_node= g.add_component(vid_metamer, label='Egreen', edge_type='/',length=0.,po=args['Epo'], diam=args['Ed'] )
                assert edge_type == '/'
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
        plant = dp['plant'][i]
        axe = dp['axe'][i]
        num_metamer = dp['numphy'][i]

        #plant, axe, num_metamer = [int(convert(d.get(x),undef=None)) for x in topology]
        print 'plant: %d, axe:%d, nb_metamers: %d'%(plant, axe, num_metamer)
        # Add new plant
        if plant != prev_plant:
            label = 'plant'+str(plant)
            vid_plant = g.add_component(g.root, label=label)
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
                vid_axe = g.add_child(vid_axe, edge_type=edge_type, label=label)
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
            vid_metamer, vid_axe =  g.add_child_and_complex(metamers[axe-1], complex=vid_axe, edge_type='+',label=label, **args)
            vid_node = nodes[axe-1]
        else:
            edge_type = '<'
            vid_metamer = g.add_child(vid_metamer, label=label, edge_type='<',**args)

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
                    vid_node= g.add_component(vid_metamer, label='Esen', edge_type='/',
                        length=args['Ev'],po=args['Epos'], diam=args['Ed'] )
                    assert edge_type == '/'
                else:
                    vid_node, vid_metamer= g.add_child_and_complex(vid_node, 
                        complex=vid_metamer, label='Esen', edge_type=edge_type,
                        length=args['Ev'],po=args['Epos'], diam=args['Ed'] )
            else:
                if vid_node == -1:
                    vid_node= g.add_component(vid_metamer, label='Egreen', edge_type='/',
                        length=args['Ev']-args['Esen'],po=args['Epo'], diam=args['Ed'] )
                    assert edge_type == '/'
                else:
                    vid_node, vid_metamer= g.add_child_and_complex(vid_node, 
                            complex=vid_metamer, label='Egreen', edge_type=edge_type,
                            length=args['Ev']-args['Esen'],po=args['Epo'], diam=args['Ed'] )
                vid_node = g.add_child(vid_node, label='Esen', edge_type='<',
                    length=args['Esen'],po=args['Epos'], diam=args['Ed'])
        else:
            if vid_node == -1:
                vid_node= g.add_component(vid_metamer, label='Egreen', edge_type='/',length=0.,po=args['Epo'], diam=args['Ed'] )
                assert edge_type == '/'
            else:
                vid_node, vid_metamer= g.add_child_and_complex(vid_node, complex=vid_metamer, label='Egreen', edge_type=edge_type,length=0.,po=args['Epo'], diam=args['Ed'] )
        # Sheath
        if args['Gv'] > 0.:
            if args['Gsen'] > args['Gv']:
                vid_node= g.add_child(vid_node, label='Gsen', edge_type='<',
                    length=args['Gv'],po=args['Gpos'], diam=args['Gd'] )
            else:
                vid_node = g.add_child(vid_node, label='Ggreen', edge_type='<',
                    length=args['Gv']-args['Gsen'],po=args['Gpo'], diam=args['Gd'] )
                vid_node = g.add_child(vid_node, label='Gsen', edge_type='<',
                    length=args['Gsen'],po=args['Gpos'], diam=args['Gd'])

        if axe == 0:
            nodes.append(vid_node)

        # Laminae
        if args['Lv'] > 0.:
            if args['Lsen'] > args['Lv']:
                l_node= g.add_child(vid_node, label='Lsen', edge_type='+',length=args['Lv'],
                    po=args['Lpos'], Ll=args['Ll'], Lw=args['Lw'], 
                    LcType=args['LcType'], LcIndex=args['LcIndex'], Linc=args['Linc'], 
                    Laz=args['Laz'], srb=0, srt=1)
            else:
                l_node = g.add_child(vid_node, label='Lgreen', edge_type='+',
                    length=args['Lv']-args['Lsen'],po=args['Lpo'],
                    Ll=args['Ll'], Lw= args['Lw'], LcType=args['LcType'], 
                    LcIndex=args['LcIndex'], Linc=args['Linc'], Laz=args['Laz'],
                    srb=0, srt=1-(args['Lsen']/args['Lv']))
                l_node = g.add_child(l_node, label='Lsen', edge_type='<',
                    length=args['Lsen'],po=args['Lpos'],
                    Ll=args['Ll'], Lw= args['Lw'], LcType=args['LcType'], 
                    LcIndex=args['LcIndex'], Linc=args['Linc'], Laz=args['Laz'],
                    srt=1, srb=1-(args['Lsen']/args['Lv']))

        
        prev_plant = plant
        prev_axe = axe
        prev_metamer = num_metamer

    return fat_mtg(g)


def _check_adel_parameters( params):
    return True

