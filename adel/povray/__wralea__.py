
# This file has been generated at Tue Feb 12 14:41:01 2008

from openalea.core import *


__name__ = 'alinea.adel.povray'
__alias__ = []

__editable__ = True
__version__ = '0.9.1'
__description__ = 'Use of povray to compute coverage rate'
__license__ = 'CECILL'
__authors__ = 'C. Pradal'
__url__ = ''
__institutes__ = 'CIRAD, INRIA'
 

__all__ = []

color_list=[(0,0,0),
            (255,0,0),
            (0,255,0),
            (0,0,255),
            (255,255,0),
            (0,255,255),
            (255,0,255),
            (128,255,0),
            (0,128,255),
            (255,0,128),
            (0,255,128),
            (128,0,255),
            (255,128,0),
            (128,128,255),
            (255,128,128),
            (128,255,128),
            (255,255,255)
            ]


povray = Factory(name='povray', 
                category='image', 
                nodemodule='povray',
                nodeclass='povray',
                inputs=[dict(name='scene', desc='PlantGL 3D scene'),
                dict(name='pov_file', interface=IStr, value='./scene.pov'),
                dict(name='camera_distance', interface=IFloat(step=1), value=200),
                dict(name='fov', interface=IInt, value=45),
                dict(name='width', interface=IInt, value=320),
                dict(name='height', interface=IInt, value=280),
                dict(name='domain', interface=ITuple, value=((0,0),(1,1))),
                dict(name='azimuth', interface=IInt, value=0),
                dict(name='zenith', interface=IInt, value=0),
                dict(name='camera_type', 
                     interface=IEnumStr(enum=['perspective', 'orthographic', 'fisheye']), 
                     value='perspective'),
                dict(name='soil', interface=IBool, value=False),
                dict(name='command', interface=IStr, value='povray'),

                ],
                outputs=[{'interface': IFileStr, 'name': 'povray image'},
                {'interface': IFileStr, 'name': 'stand box image'}],
                )

__all__.append('povray')

color_item = Factory(name='col_item', 
                category='color', 
                nodemodule='povray',
                nodeclass='col_item',
                inputs=[dict(interface='IInt', name='color index', value=None, desc='color index. If None, return a function'),
                        dict(interface='ISequence', name='color list', value=color_list),
                        ],
                outputs=[{'interface': IRGBColor, 'name': 'color list'}],
                )

__all__.append('color_item')


domain3D = Factory(name='domain3D', 
                category='image', 
                nodemodule='povray',
                nodeclass='domain3D',
                inputs=[dict(interface='ITuple', name='domain2D', value=()),
                        dict(interface=None, name='scene', value=None),
                        ],
                outputs=[{'interface': ITuple, 'name': 'domain3D'}],
                )

__all__.append('domain3D')


stand_box = Factory(name='stand_box', 
                category='image', 
                nodemodule='povray',
                nodeclass='stand_box',
                inputs=[dict(interface='ITuple', name='domain', value=((0, 0, 0), (1, 1, 1)))],
                outputs=[{'interface': None, 'name': 'stand_box'}],
                )

__all__.append('stand_box')


