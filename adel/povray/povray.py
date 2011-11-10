import os
from openalea.core.path import path

def povray(pov_file, height=320, width=280, image_name=''):
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

def povray2(scene, pov_file = '', camera_distance=1., fov=45., width=320, height=280, domain = ((0.,0.),(1.,1.)), azimuth=0., zenith= 0., image_name=''):
    '''    
    Compute povray file based on a scene.
    :Parameters:
        - scene: a plantgl scene
        - camera_distance: distance from the position of the camera to the look_at point
        - angle : angle corresponding to the width of the image
        - width: width of the final image in pixel
        - height: height of the final image in pixel
        - domain: scene pattern used by caribu
        - azimuth: angle in degree around the vertical axis. az=0. is equialent to have the width of the image align to X direction of the scene
        - zenith: angle between the view direction and the vertical
    '''

    # color elements based on the optical properties
    
    t=Tesselator()
    mesh_fn = pov_file[0:-4]+'_mesh.pov'
    pov = PovFilePrinter(mesh_fn, t)
    scene.apply(pov)

    # Write the camera in the main pov file
    #pov_camera = pov_header(...)
    pov_camera = """
camera {
    perspective
    location <0,0,{z}>
    direction <0,0,-1>
    right <36/24,0,0>
    look_at <0,0,0 >
    angle {fov}

    rotate <0,0,{azimuth}>
    rotate <{zenith},0,0>
    translate <{x},{y},0>
}

    """

    pov_camera.format(z=camera_distance, fov=fov, azimuth=azimuth, zenith=zenith, 
        x=0.5*(domain[0][0]+domain[1][0]), y=0.5*(domain[0][1]+domain[1][1]))
    
    f = path(pov_file)
    namebase = f.namebase
    ext = f.ext
    dirname = f.dirname()

    povf = f.open(mode='w')
    pov_header(povf, pov_camera)
    povf.write('#include "{mesh_fn}"\n'.format(mesh_fn=mesh_fn))
    povf.close()

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

def pov_header(f, pov_camera, background = (0,0,0), light = (0,0,10)):
    """ Write the header of  povray file. """
    f.write("/"*80+"\n")
    f.write("#include \"colors.inc\""+"\n\n")
    f.write(pov_camera+"\n")
    x,y,z=light
    f.write("light_source {  <{x},{y},{z}> color rgb <1,1,1>}\n\n".format(x=x,y=y,z=z))
    r,g,b = background
    f.write("background { color rgb < {r},{g},{b}> }\n\n".format(r=r,g=g,b=b))

