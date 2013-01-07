""" Methods for mtg interpretation with turtle """

from math import degrees, radians, pi, cos, sin
import openalea.plantgl.all as pgl
from openalea.mtg import turtle, MTG
from openalea.mtg.algo import union


# Meshing function for StemElements

def slim_cylinder(length, radius_base, radius_top):
    " Try to construct a cylinder with a low number of triangles. "

    rb, rt = radius_base, radius_top
    a1, a2, a3 = 0, 2*pi/3., 4*pi/3.
    r = rb
    p1 = (r*cos(a1), r*sin(a1),0)
    p2 = (r*cos(a2), r*sin(a2),0)
    p3 = (r*cos(a3), r*sin(a3),0)
    r = rt
    q1 = (r*cos(a1+pi), r*sin(a1+pi),length)
    q2 = (r*cos(a2+pi), r*sin(a2+pi),length)
    q3 = (r*cos(a3+pi), r*sin(a3+pi),length)
    set = pgl.TriangleSet([p1, p2, p3, q1, q2, q3],
                      [(2,1,0), (3,4,5), (0,5,4), (0,4,2), (2,4,3), (3,1,2), (1,3,5), (5,0,1)])
    return set
    
def StemElement_mesh(length, diameter_base, diameter_top, classic = False):
    """ Compute mesh for a stem element
        - classic indicates
    """
    if classic:
        solid = True
        slices = 3
        diameter = diameter_base
        stem = pgl.Tapered(diameter_base/2., diameter_top/2., pgl.Cylinder(1., length , solid, slices))
        tessel = pgl.Tesselator()
        stem.apply(tessel)
        mesh = tessel.triangulation
    else:
        mesh = slim_cylinder(length, diameter_base /2., diameter_top /2.)

    return mesh

#meshing function for leaf elements 
    
def LeafElement_mesh(shape, L_shape, Lw_shape, length, s_base, s_top):
    """ Compute mesh for a leaf element.
        - shape is a x,y,s,r tuple descriibing leaf shape
        - L_shape is the length of the scaled shape
        - Lw_shape is the width of the scaled shape
        - length is the total visible length to be meshed
        - s_base and s_top are relative proportion (on length) of the element to represent
    """
    
    from alinea.adel.fitting import mesh4, plantgl_shape
    
    leaf_mesh = mesh4(shape, L_shape, length, s_base, s_top, Lw_shape)
    if leaf_mesh:
        pts, ind = leaf_mesh
        if len(ind) < 2:
            mesh = None
        else:
            mesh = plantgl_shape(pts, ind)
    else:
        mesh = None

    return mesh

    
def incline_leaf(shape, inclin, relative_angle = True):
    """ transform a xysr tuple representing leaf shape to get a given angle at leaf base.
     - angle the desired angle (deg)
     - if relative_angle == True, angle is interpreted as a multiplier to original shape angle
     """   
    Linc = inclin
    x, y = shape[0], shape[1]
    init_angle = pgl.angle((x[1]-x[0], y[1]-y[0]),(0,1))

    if relative_angle:
        angle = Linc * init_angle
        angle = min(pi, angle)
    else:
        angle = radians(Linc)
    
    rotation_angle = init_angle - angle

    # rotation of the midrib
    cos_a = cos(rotation_angle); sin_a = sin(rotation_angle)

    x1 = x[0] + cos_a*x - sin_a*y
    y1 = y[0] + sin_a*x + cos_a*y
    leaf= x1, y1, shape[2], shape[3]
    
    return leaf
    
def compute_element(element_node, classic=False): 
    """ compute geometry of Adel base elements (LeafElement and StemElement) 
    element_node should be a mtg node proxy"""
    n = element_node
    geom = None
    
    if n.label.startswith('Leaf'): #leaf element
        blade = n.complex()
        if blade.shape_xysr:
            if blade.inclination > 0:
                shape = incline_leaf(blade.shape_xysr, blade.inclination)
            else: 
                shape = blade.shape_xysr
            # x-> -x as turtlle expect x-y+
            leaf = (-shape[0],)+shape[1:]
            geom = LeafElement_mesh(leaf, blade.shape_mature_length, blade.shape_max_width, 
                                n.length, n.srb, n.srt)   
    elif n.label.startswith('Stem'): #stem element
        stem = n.complex()
        #diameter_base = stem.parent().diameter if (stem.parent() and stem.parent().diameter > 0.) else stem.diameter
        #diameter_top = n.diam
        diameter_base = stem.diameter
        diameter_top = stem.diameter
        geom = StemElement_mesh(n.length, diameter_base, diameter_top, classic)
        
    return geom
        
def _transform(turtle, mesh):
        x = turtle.getUp()
        z = turtle.getHeading()

        bo = pgl.BaseOrientation(x, z^x)
        matrix = pgl.Transform4(bo.getMatrix())
        matrix.translate(turtle.getPosition())
        #print 'Position ', turtle.getPosition()
        mesh = mesh.transform(matrix)
        return mesh

def adel_visitor(g, v, turtle):
    """ Performs geometric interpretation of mtg nodes
    """
    # 1. retrieve the node
    n = g.node(v)

    # Go to plant position if first plant element
    if n.parent() is None:#this is a new plant base
        p = n.complex_at_scale(scale=1)
        if 'position' in p.properties():
            #print 'moving to ', p.position
            turtle.move(p.position)
        else:
            turtle.move(0,0,0)
        turtle.setHead(0,0,1,1,0,0)
        if 'azimuth' in p.properties():
            turtle.rollL(p.azimuth)
    #hypothesis that inclin is to be applied at the base of the visible elements
    if n.offset > 0:
        turtle.f(n.offset)
    #incline turtle at the base of stems,
    if n.label.startswith('Stem'):
        inclin = float(n.inclination) if n.inclination else 0.
        azim = float(n.azimuth) if n.azimuth else 0.
        if inclin:
            # incline along curent azimuth for ramification (tiller bases)
            if n.edge_type() == '+':
                #print ' + node ', n._vid, 'up before inclin ', turtle.getUp()
                turtle.down(inclin)
            # if not incline towardss vertical
            else:
                up = turtle.getUp()
                turtle.rollToVert()
                angle = degrees(pgl.angle(up,turtle.getUp()))
                turtle.down(inclin)
                #replace turtle in original azimuth plane
                #print 'angle ', angle
                turtle.rollL(-angle)
        if azim:
            #print 'node', n._vid, 'azim ', azim
            turtle.rollL(azim)

       
    #if n.label.startswith('Leaf'):    
    
    # update geometry of elements
    if n.length > 0:
        mesh = compute_element(n)
        if mesh:#To DO : reset to None if calculated so ?
            n.geometry = _transform(turtle, mesh)
    # 3. Update the turtle
    turtle.setId(v)
    if n.label.startswith('Stem'):
        if n.length > 0:
            turtle.f(n.length)
   




    
        
def mtg_interpreter(g, visitor = adel_visitor):
    ''' Compute/update the geometry on each node of the MTG using Turtle geometry. '''
#BUG : sub_mtg mange le vertex plant => on perd la plante !
    #plants = g.component_roots_at_scale(g.root, scale=1)
    #nplants = g.nb_vertices(scale=1)
    #gt = MTG()
    
    #for plant in plants:
    #   gplant = g.sub_mtg(plant)
    scene = turtle.TurtleFrame(g,visitor=visitor, gc=False, all_roots=True)
    #   gt = union(gplant,gt)
       
    return g
       
def plot3d(g, 
               leaf_material = None,
               stem_material = None,
               soil_material = None,
               colors = None):
    """
    Returns a plantgl scene from an mtg.
    """
    
    Material = pgl.Material
    Color3 = pgl.Color3
    Shape = pgl.Shape
    Scene = pgl.Scene
    
    if colors is None:
        if leaf_material is None:
            leaf_material = Material(Color3(0,180,0))
        if stem_material is None:
            stem_material = Material(Color3(0,130,0))
        if soil_material is None:
            soil_material = Material(Color3(170, 85,0))
        colors = g.property('color')
        
    geometries = g.property('geometry')
    greeness = g.property('is_green')
    labels = g.property('label')
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
        is_green = greeness.get(vid)
        if colors:
            shape = Shape(mesh, Material(Color3(* colors.get(vid, [0,0,0]) )))
        elif not greeness:
            if not shape:
                shape = Shape(mesh)
        elif label.startswith('Stem') and is_green:
            shape = Shape(mesh, stem_material)
        elif label.startswith('Leaf') and is_green:
            shape = Shape(mesh, leaf_material)
        elif not is_green:
            shape = Shape(mesh, soil_material)
        shape.id = vid
        scene.add(shape)

    for vid, mesh in geometries.iteritems():
        geom2shape(vid, mesh, scene)
    return scene

        