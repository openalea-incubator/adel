
# This file has been generated at Fri Jun 25 13:19:57 2010

from openalea.core import *


__name__ = 'alinea.prospect'

__editable__ = True
__description__ = ''
__license__ = None
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = None
__authors__ = ''
__institutes__ = None
__icon__ = 'leaf.png'


__all__ = ['prospect_prospect_5', '_163038604']



prospect_prospect_5 = Factory(name='prospect',
                description='prospect compute the radiative optical leaf properties.',
                category='Unclassified',
                nodemodule='prospect',
                nodeclass='prospect_5',
                inputs=[{'interface': IFloat, 'name': 'N', 'value': 1.0}, {'interface': IFloat, 'name': 'Cab', 'value': 0.0}, {'interface': IFloat, 'name': 'Car', 'value': 0.0040000000000000001}, {'interface': IFloat, 'name': 'Cw', 'value': 0.0040000000000000001}, {'interface': IFloat, 'name': 'Cm', 'value': 0.0019}],
                outputs=({'name': 'wavelenght'}, {'name': 'reflectance'}, {'name': 'transmittance'}),
                widgetmodule=None,
                widgetclass=None,
               )




_163038604 = CompositeNodeFactory(name='Demo Corn',
                             description='',
                             category='Unclassified',
                             doc='',
                             inputs=[],
                             outputs=[],
                             elt_factory={  2: ('alinea.prospect', 'prospect'),
   3: ('openalea.pylab.nodes', 'PyLabScatter'),
   10: ('openalea.data structure', 'float scy')},
                             elt_connections={  148622796: (2, 0, 3, 0), 148622820: (10, 0, 2, 4), 148622844: (2, 1, 3, 1)},
                             elt_data={  2: {  'block': False,
         'caption': 'prospect',
         'factory': '<openalea.core.node.NodeFactory object at 0x9b7c02c> : "prospect"',
         'hide': True,
         'id': 2,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': -164.0,
         'posy': -18.0,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   3: {  'block': False,
         'caption': 'PyLabScatter',
         'factory': '<openalea.core.node.NodeFactory object at 0x9d797cc> : "PyLabScatter"',
         'hide': True,
         'id': 3,
         'lazy': False,
         'port_hide_changed': set(),
         'posx': -250.83771720070143,
         'posy': 48.675753228120527,
         'priority': 0,
         'use_user_color': False,
         'user_application': False,
         'user_color': None},
   10: {  'block': False,
          'caption': '3.7e-03',
          'factory': '<openalea.core.node.NodeFactory object at 0xa8e9c0c> : "float scy"',
          'hide': True,
          'id': 10,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': -148.0,
          'posy': -108.0,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'hide': True,
                'id': 0,
                'lazy': True,
                'port_hide_changed': set(),
                'posx': 0,
                'posy': 0,
                'priority': 0,
                'use_user_color': True,
                'user_application': None,
                'user_color': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'hide': True,
                 'id': 1,
                 'lazy': True,
                 'port_hide_changed': set(),
                 'posx': 0,
                 'posy': 0,
                 'priority': 0,
                 'use_user_color': True,
                 'user_application': None,
                 'user_color': None}},
                             elt_value={  2: [(0, '1.5'), (1, '58.0'), (2, '20.0'), (3, '0.01')],
   3: [  (2, '20'),
         (3, "'blue'"),
         (4, "'circle'"),
         (5, 'None'),
         (6, '0.5'),
         (7, 'True'),
         (8, 'True'),
         (9, '1'),
         (10, "'Wave Length'"),
         (11, "'Reflectance'"),
         (12, "''"),
         (13, '1'),
         (14, 'True'),
         (15, 'False'),
         (16, '{}'),
         (  17,
            "{'xmin': None, 'ymin': None, 'type': 'normal', 'ymax': None, 'xmax': None}")],
   10: [(0, "'0.003662'")],
   '__in__': [],
   '__out__': []},
                             elt_ad_hoc={  2: {'position': [-164.0, -18.0], 'userColor': None, 'useUserColor': False},
   3: {'position': [-250.83771720070143, 48.675753228120527], 'userColor': None, 'useUserColor': False},
   4: {  'position': [-265.0, 66.0], 'useUserColor': False, 'userColor': None},
   5: {  'position': [-208.0, -57.0], 'useUserColor': False, 'userColor': None},
   6: {  'position': [-182.0, 67.0], 'useUserColor': False, 'userColor': None},
   7: {  'position': [-105.0, 74.0], 'useUserColor': False, 'userColor': None},
   8: {  'position': [-57.0, -52.0], 'useUserColor': False, 'userColor': None},
   9: {  'position': [4.0, -26.0], 'useUserColor': False, 'userColor': None},
   10: {'position': [-148.0, -108.0], 'userColor': None, 'useUserColor': False},
   '__in__': {'position': [0, 0], 'userColor': None, 'useUserColor': True},
   '__out__': {'position': [0, 0], 'userColor': None, 'useUserColor': True}},
                             lazy=True,
                             )




