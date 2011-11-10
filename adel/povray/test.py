"""
from povray import *

scene = Scene()
scene.add(Sphere(radius=0.1))

povray(scene)

"""

pov_camera = """
camera {{
    perspective
    location <0,0,{0:f}>
    direction <0,0,-1>
    right <36/24,0,0>
    look_at <0,0,0 >
    angle {1:f}

    rotate <0,0,{2:f}>
    rotate <{3:f},0,0>
    translate <{4:f},{5:f},0>
}}
"""

#pov_camera = pov_camera.format(tz=1., fov=1., azimuth=1.,  zenith=1., tx=1., ty=1.)
pov_camera = pov_camera.format(1., 1., 1.,  1., 1., 1.)
