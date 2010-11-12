
# This file has been generated at Thu Oct 28 23:59:19 2010

from openalea.core import *


__name__ = 'alinea.adel.geometry'

__editable__ = True
__description__ = 'Geometric features for plant reconstruction and scene creation'
__license__ = 'CECILL'
__url__ = ''
__alias__ = ['adel.geometry']
__version__ = '0.0.1'
__authors__ = 'C. Pradal'
__institutes__ = 'CIRAD'
__icon__ = ''


__all__ = ['geometry_symbols', 'geometry_mtg_turtle', 'geometry_leaf_to_mesh', 'setGeometry_setGeometry', 'geometry_leaf_element']



geometry_symbols = Factory(name='symbols',
                authors='C. Pradal (wralea authors)',
                description='Build symbols like Leaf and stem for Turtle interpretation',
                category='graphic',
                nodemodule='geometry',
                nodeclass='symbols',
                inputs=[{'interface': IDict, 'name': 'leaves', 'value': {}, 'desc': 'database'}, {'interface': IInt, 'name': 'seed', 'value': None, 'desc': 'random seed'}],
                outputs=[{'interface': IDict, 'name': 'symbols', 'desc': 'A set of symbols'}],
                widgetmodule=None,
                widgetclass=None,
               )




geometry_mtg_turtle = Factory(name='MTG Interpreter',
                authors='C. Pradal (wralea authors)',
                description='Geometric embedding of the MTG',
                category='geometry',
                nodemodule='geometry',
                nodeclass='mtg_turtle',
                inputs=[{'interface': None, 'name': 'g', 'value': None, 'desc': ''}, {'interface': None, 'name': 'symbols', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'g', 'desc': 'MTG with geometry'}],
                widgetmodule=None,
                widgetclass=None,
               )




geometry_leaf_to_mesh = Factory(name='leaf to mesh',
                authors='C. Pradal (wralea authors)',
                description='convert a midrib data to a mesh',
                category='graphic',
                nodemodule='geometry',
                nodeclass='leaf_to_mesh',
                inputs=[{'interface': ISequence, 'name': 'leaf', 'value': [], 'desc': ''}, {'interface': IFloat, 'name': 'length_max', 'value': 1.0, 'desc': ''}, {'interface': IFloat, 'name': 'length', 'value': 1.0, 'desc': ''}, {'interface': IFloat, 'name': 'radius', 'value': 1.0, 'desc': ''}],
                outputs=[{'interface': IInterface, 'name': 'mesh', 'desc': 'A plantgl mesh'}],
                widgetmodule=None,
                widgetclass=None,
               )




setGeometry_setGeometry = Factory(name='setGeometry',
                authors='C. Pradal (wralea authors)',
                description='Control panel for geometry settings',
                category='data i/o',
                nodemodule='setGeometry',
                nodeclass='setGeometry',
                inputs=[{'interface': IInt, 'name': 'Poygons per leaf', 'value': None, 'desc': ''}, {'interface': ITextStr, 'name': 'axe azimth', 'value': '#This function should return axe azimuth as a function of axe number a (a =0 for bm)\nfunction(a) {\nifelse(a==0,0,75)\n}', 'desc': ''}, {'interface': ITextStr, 'name': 'axe inclination', 'value': None, 'desc': ''}, {'interface': ITextStr, 'name': 'axe distance to row', 'value': None, 'desc': ''}, {'interface': ITextStr, 'name': 'leaf azimuth', 'value': None, 'desc': ''}, {'interface': ITextStr, 'name': 'leaf index', 'value': None, 'desc': ''}],
                outputs=[{'interface': IInt, 'name': 'polygons', 'desc': ''}, {'interface': None, 'name': 'geomAxe', 'desc': 'R list foraxe geometry'}, {'interface': None, 'name': 'geomLeaf', 'desc': 'R list for leaf geometry'}],
                widgetmodule=None,
                widgetclass=None,
               )




geometry_leaf_element = Factory(name='leaf element',
                authors='C. Pradal (wralea authors)',
                description='convert a midrib data to a mesh',
                category='graphic',
                nodemodule='geometry',
                nodeclass='leaf_element',
                inputs=[{'interface': ISequence, 'name': 'leaf', 'value': [], 'desc': ''}, {'interface': IFloat, 'name': 'length_max', 'value': 1.0, 'desc': ''}, {'interface': IFloat, 'name': 'length', 'value': 1.0, 'desc': ''}, {'interface': IFloat, 'name': 's_min', 'value': 0.0, 'desc': ''}, {'interface': IFloat, 'name': 's_max', 'value': 1.0, 'desc': ''}, {'interface': IFloat, 'name': 'radius', 'value': 1.0, 'desc': ''}],
                outputs=[{'interface': IInterface, 'name': 'mesh', 'desc': 'A plantgl mesh'}],
                widgetmodule=None,
                widgetclass=None,
               )




