from math import *


def lightString(radiance, zenith_angle, azimuth_angle):
    '''    compute the directinal vactor from zenith and azimuth angles
    '''
    theta = zenith_angle / 180. * pi
    phi = azimuth_angle / 180. * pi
    vector = (radiance, sin(theta) * cos(phi),sin(theta) * sin(phi),  - cos(theta)) 
    return  ' '.join(map(str,vector))
