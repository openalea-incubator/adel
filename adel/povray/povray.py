import os
import platform
from openalea.core.path import path
from openalea.plantgl.all import *

windows = platform.system() is 'Windows'

def povray2(pov_file, height=320, width=280, image_name=''):
    '''    
    '''
    image = None; 
    # write the node code here.

    f = path(pov_file)
    namebase = f.namebase
    ext = f.ext
    dirname = f.dirname()

    old_dir = os.path.abspath(os.getcwd())
    os.chdir(dirname)

    if not image_name:
        image_name = dirname / namebase +'.png'

    image_name = path(image_name)
    if image_name.exists():
        image_name.remove()

    cmdline = 'povray +I%s +H%d +W%d'%(pov_file, height, width)
    os.system(cmdline)
    # return outputs
    os.chdir(old_dir)

    return image_name,

def povray(scene, 
           pov_file = './scene.pov', 
           camera_distance=1., fov=45., width=320, height=280, 
           domain = ((0.,0.),(1.,1.)), 
           azimuth=0., zenith= 0., camera_type = 'perspective', 
           soil=False, 
           povray_cmd='povray'):
    """    
    Compute povray files based both on a scene and its stand box.

    :Parameters:
        - scene: a plantgl scene
        - camera distance: distance from the position of the camera to the look_at point
        - angle : angle corresponding to the width of the image
        - width: width of the final image in pixel
        - height: height of the final image in pixel
        - domain: scene pattern used by caribu
        - azimuth: angle in degree around the vertical axis. az=0. is equialent to have the width of the image align to X direction of the scene
        - zenith: angle between the view direction and the vertical
        - camera type: perspective, orthographic or fisheye
        - soil: add a soil to the scene
        - povray cmd : the path of the povray exe
    """

    # color elements based on the optical properties
    
    t=Tesselator()
    mesh_fn = pov_file[0:-4]+'_mesh.pov'
    pov = PovFilePrinter(mesh_fn, t)
    scene.apply(pov)
    
    d3D = domain3D(domain, scene)
    
    scene_box = Scene()
    scene_box.add(stand_box(d3D))
    
    mesh_fn_box = pov_file[0:-4]+'_mesh_box.pov'
    pov_box = PovFilePrinter(mesh_fn_box, t)
    scene_box.apply(pov_box)
    
    # Write the camera in the main pov file
    #pov_camera = pov_header(...)
    pov_camera = """
camera {{
    {camera}
    location <0,0,{tz:.2f}>
    direction <0,0,-1>
    right <{width:d}/{height:d},0,0>
    look_at <0,0,0 >
    angle {fov:n}

    rotate <{zenith:.2f},0,0>
    rotate <0,0,{azimuth:.2f}>
    translate <{tx:.2f},{ty:.2f},0>
}}

    """

    x=0.5*(domain[0][0]+domain[1][0])
    y=0.5*(domain[0][1]+domain[1][1])
    pov_camera = pov_camera.format(camera=camera_type, tz=camera_distance, fov=fov, azimuth=azimuth, 
                      zenith=zenith, tx=x, ty=y, width=width, height=height)
    
    #print pov_camera

    f = path(pov_file)
    namebase = f.namebase
    ext = f.ext
    dirname = f.dirname()
    f_box = dirname / namebase + '_box' + ext  

    povf = f.open(mode='w')
    pov_header(povf, pov_camera, domain=domain, soil=soil)
    povf.write('#include "{mesh_fn}"\n'.format(mesh_fn=mesh_fn))
    povf.close()
    
    povf_box = f_box.open(mode='w')
    pov_header(povf_box, pov_camera, domain=domain,light = (x,y,camera_distance))
    povf_box.write('#include "{mesh_fn}"\n'.format(mesh_fn=mesh_fn_box))
    povf_box.close()

    old_dir = os.path.abspath(os.getcwd())
    os.chdir(dirname)

    #if windows:
    #        image_name = dirname / namebase + '.bmp'
    #        image_name_box = dirname / namebase + '_box' + '.bmp'
    #else:
    #    image_name = dirname / namebase + '.png'
    #    image_name_box = dirname / namebase + '_box' + '.png'
    
    image_name = dirname / namebase + '.png'
    image_name_box = dirname / namebase + '_box' + '.png'

    image_name_box = path(image_name_box)
    image_name = path(image_name)
    if image_name.exists():
        image_name.remove()
    if image_name_box.exists():
        image_name_box.remove()
    
    pov_file_box = pov_file[0:-4] +'_box.pov'
    if windows:
        cmdline = 'pvengine +I%s +H%d +W%d -d /exit'%(pov_file, height, width)
        cmdline_box = 'pvengine +I%s +H%d +W%d -d /exit'%(pov_file_box, height, width)
    else:
        cmdline = '%s +I%s +H%d +W%d'%(povray_cmd, pov_file, height, width)
        cmdline_box = '%s +I%s +H%d +W%d'%(povray_cmd, pov_file_box, height, width)
    
    os.system(cmdline)
    os.system(cmdline_box)

    if not image_name.exists():
        image_name = dirname / namebase + '.bmp'
        image_name_box = dirname / namebase + '_box' + '.bmp'

    if not os.path.isfile(image_name):
        print 'Error: Image not created. Check that povray is installed and set in the path.'
    if not os.path.isfile(image_name_box):
        print 'Error: Image box not created. Check that povray is installed and set in the path.'
    # return outputs
    os.chdir(old_dir)

    return image_name, image_name_box

    
def pov_header(f, pov_camera, background = (0,0,0), light = (10,0,10), domain = ((0.,0.),(1.,1.)), soil=False,light_color=(1,1,1)):
    """ Write the header of  povray file. """
    f.write("/"*80+"\n")
    f.write("#include \"colors.inc\""+"\n\n")
    f.write(pov_camera+"\n")
    x,y,z=light
    f.write("light_source {{  <{tx},{ty},{tz}> color rgb <{r},{g},{b}>}}\n\n".format(tx=x,ty=y,tz=z,r=light_color[0],g=light_color[1],b=light_color[2]))
    r,g,b = background
    f.write("background {{ color rgb < {cr},{cg},{cb}> }}\n\n".format(cr=r,cg=g,cb=b))
    
    if soil:
        (x1,y1),(x2,y2) = domain
        s = """
#declare T1=
    texture {{
      pigment {{ color red 1.0 green 0.75 blue 0.33 }}
      normal {{bumps 0.4  scale 0.007}} 
}}

#declare T2=
 texture{{
   pigment{{color red 1.0 green 0.75 blue 0.33}}
   normal {{granite 0.6 scale 0.1 }}
   finish {{phong 0.8 phong_size 200}} 
}}

box {{ <{x1}, {y1},  -0.1>,  
    < {x2}, {y2}, 0>   
    texture {{T2}} }}
"""
        text = s.format(x1=x1,y1=y1,x2=x2,y2=y2)
        print text
        f.write(text)
        

def domain3D(domain2D, scene):
    t=Tesselator()
    bbc = BBoxComputer(t)
    bbc.process(scene)
    bbox = bbc.result

    z_base = bbox.getZMin()
    z_top = bbox.getZMax()
    domain3D = (domain2D[0] + (z_base,), domain2D[1] + (z_top,))
    return domain3D
        
        
def stand_box(domain):
    '''
    
    domain: 3D bounding box of the stand
    '''
    # list of points
    z_base = domain[0][2]
    z_top = domain[1][2]
    sides_points = [domain[0], # coordinates of bottom right corner
                  (domain[1][0], domain[0][1], z_base),
                  (domain[1][0], domain[1][1], z_base), # coordinates of bottom left corner
                  (domain[0][0], domain[1][1], z_base),
                  (domain[0][0], domain[0][1], z_top), # coordinates of top right corner
                  (domain[1][0], domain[0][1], z_top),
                  (domain[1][0], domain[1][1], z_top),    # coordinates of top left corner
                  (domain[0][0], domain[1][1], z_top)]
                  
    bottom_points = [domain[0], # coordinates of bottom right corner
                  (domain[1][0], domain[0][1], z_base),
                  (domain[1][0], domain[1][1], z_base), # coordinates of bottom left corner
                  (domain[0][0], domain[1][1], z_base)]

    # list of indices to make the quads of the sides from the points
    side_indices = [(0, 1, 5, 4), #
               (1, 2, 6, 5), # indices for 
               (2, 3, 7, 6), # side faces
               (3, 0, 4, 7)] #         
     
    # list of indices to make the quads of the bottom from the points
    bottom_indices = [(0, 1, 2, 3)] # indices for bottom face

    # list of colors
    side_color = Color3(0, 0, 0)
    bottom_color = Color3(255, 255, 255)

    # construction of the geometry for the sides
    side_box = QuadSet(sides_points, side_indices)
    # construction of the geometry for the bottom
    bottom_box = QuadSet(bottom_points, bottom_indices)
                           
    # create 2 shapes: 1 with side_color, 1 with bottom_color 
    sides_shape = Shape(side_box, Material(side_color))
    bottom_shape = Shape(bottom_box, Material(bottom_color))
    
    scene = Scene()
    scene.add(sides_shape)
    scene.add(bottom_shape)
    
    return scene
    
    
color_list=[(0,0,0),
            (255,0,0),# 1 = Green lamina
            (0,255,0),# 2 = Senescent Lamina
            (0,0,255),# 3 = Green sheath
            (255,255,0),# 4 = Senescent sheath
            (0,255,255),# 5 = Green internode
            (255,0,255),# 6 = senescent Internode
            (128,255,0),# 7 = green peduncule
            (0,128,255),# 8 = senescent peduncule
            (255,0,128),# 9 = green ear
            (0,255,128),# 10 = senescent ear
            (128,0,255),# 11 = green awn
            (255,128,0),# 12 = senescent awn
            (128,128,255),#???
            (255,128,128),
            (128,255,128),
            (255,255,255)
            ]

def col_item (ind, color_list=color_list) :
    if ind is None:
        return lambda x: color_list[(x-1) % len(color_list)]
    else:
        return color_list[(ind-1) % len(color_list)],

        