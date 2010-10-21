
# This file has been generated at Fri Dec 12 09:25:36 2008

from openalea.core import *


__name__ = 'alinea.caribu.visualisation'

__editable__ = True
__description__ = ' Visualisation tools for the Caribu package '
__license__ = ''
__url__ = ''
__alias__ = ['Caribu.Visualisation']
__version__ = '0.0.3'
__authors__ = 'C. Pradal'
__institutes__ = 'CIRAD'
__icon__ = ''
 

__all__ = ['py_canview_read_can', 'gammaTrans_gammaTrans', 'ViewMapOnCan', '_161866188', 'py_canview_plot_can']



py_canview_read_can = Factory(name='Import Can File',
                description='ld a detailled description of a can file',
                category='IO, Visualisation.Caribu',
                nodemodule='py_canview',
                nodeclass='read_can',
                inputs=[{'interface': IFileStr(filter="*.can", save=False), 'name': 'can file'}],
                outputs=[{'name': 'Canestra Scene'}],
                widgetmodule=None,
                widgetclass=None,
               )




gammaTrans_gammaTrans = Factory(name='gammaTrans',
                description='return value normalised and raised at exponent gamma',
                category='Unclassified',
                nodemodule='gammaTrans',
                nodeclass='gammaTrans',
                inputs=[{'interface': None, 'name': 'values', 'value': None, 'desc': ''}, {'interface': IFloat, 'name': 'gamma', 'value': 1, 'desc': ''}],
                outputs=[{'interface': ISequence, 'name': 'res', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




ViewMapOnCan = CompositeNodeFactory(name='ViewMapOnCan', 
                             description='Plots values on a canestra scene', 
                             category='scene',
                             doc='',
                             inputs=[  {'desc': '', 'interface': IStr, 'name': 'CansceneString', 'value': None},
   {'desc': '', 'interface': ISequence, 'name': 'Values', 'value': None},
   {  'desc': '',
      'interface': IFloat,
      'name': 'gamma',
      'value': 0.20000000000000001}],
                             outputs=[{'desc': '', 'interface': IFileStr, 'name': 'can file'}],
                             elt_factory={  4: ('alinea.caribu.visualisation', 'gammaTrans'),
   5: ('openalea.color', 'colormap'),
   6: ('openalea.functional', 'map'),
   7: ('alinea.caribu.visualisation', 'Plot Can File'),
   8: ('alinea.caribu', 'WriteCan'),
   9: ('openalea.file', 'tmpnam')},
                             elt_connections={  9919108: ('__in__', 2, 4, 1),
   9919120: (7, 0, '__out__', 0),
   9919132: (5, 0, 6, 0),
   9919144: (8, 0, 7, 0),
   9919156: (6, 0, 7, 1),
   9919168: (9, 0, 8, 1),
   9919180: ('__in__', 0, 8, 0),
   9919192: ('__in__', 1, 4, 0),
   9919204: (4, 0, 6, 1)},
                             elt_data={  4: {  'block': False,
         'caption': 'gammaTrans',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([]),
         'posx': 340.27777777777777,
         'posy': -423.61111111111109,
         'priority': 0,
         'user_application': None},
   5: {  'block': False,
         'caption': 'ColorMap',
         'hide': True,
         'lazy': True,
         'minimal': False,
         'port_hide_changed': set([]),
         'posx': 270.41666666666652,
         'posy': -368.61111111111114,
         'priority': 0,
         'user_application': None},
   6: {  'block': False,
         'caption': 'map',
         'hide': True,
         'lazy': True,
         'minimal': False,
         'port_hide_changed': set([]),
         'posx': 416.11111111111103,
         'posy': -283.47222222222229,
         'priority': 0,
         'user_application': None},
   7: {  'block': False,
         'caption': 'Plot Can File',
         'hide': True,
         'lazy': True,
         'minimal': False,
         'port_hide_changed': set([]),
         'posx': 120.0,
         'posy': -227.5,
         'priority': 0,
         'user_application': None},
   8: {  'block': False,
         'caption': 'WriteCan',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([]),
         'posx': 95.0,
         'posy': -308.75,
         'priority': 0,
         'user_application': None},
   9: {  'block': False,
         'caption': 'tmpnam',
         'hide': True,
         'lazy': False,
         'minimal': False,
         'port_hide_changed': set([]),
         'posx': 155.0,
         'posy': -372.5,
         'priority': 0,
         'user_application': None},
   '__in__': {  'caption': 'In',
                'hide': True,
                'lazy': True,
                'minimal': False,
                'port_hide_changed': set([]),
                'posx': 230.0,
                'posy': -540.0,
                'priority': 0},
   '__out__': {  'caption': 'Out',
                 'hide': True,
                 'lazy': True,
                 'minimal': False,
                 'port_hide_changed': set([]),
                 'posx': 145.0,
                 'posy': -137.5,
                 'priority': 0}},
                             elt_value={  4: [],
   5: [(0, 'None'), (1, '0'), (2, '1'), (3, '250.0'), (4, '20.0')],
   6: [],
   7: [],
   8: [],
   9: [],
   '__in__': [],
   '__out__': []},
                             lazy=True,
                             )




_161866188 = CompositeNodeFactory(name='Plot CaribuScene', 
                             description='3D plot of Caribuscene', 
                             category='visualisation',
                             doc='',
                             inputs=[{'desc': '', 'interface': IStr, 'name': 'CansceneString', 'value': None}],
                             outputs=[{'desc': '', 'interface': IFileStr, 'name': 'can file'}],
                             elt_factory={  7: ('alinea.caribu.visualisation', 'Plot Can File'),
   8: ('alinea.caribu', 'WriteCan'),
   9: ('openalea.file', 'tmpnam')},
                             elt_connections={  9919168: (9, 0, 8, 1),
   9919180: (8, 0, 7, 0),
   9919192: ('__in__', 0, 8, 0),
   9919204: (7, 0, '__out__', 0)},
                             elt_data={  7: {  'block': False,
         'caption': 'Plot Can File',
         'hide': True,
         'lazy': True,
         'minimal': False,
         'port_hide_changed': set([]),
         'posx': 138.75,
         'posy': -313.75,
         'priority': 0,
         'user_application': None},
   8: {  'block': False,
         'caption': 'WriteCan',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([]),
         'posx': 130.0,
         'posy': -396.25,
         'priority': 0,
         'user_application': None},
   9: {  'block': False,
         'caption': 'tmpnam',
         'hide': True,
         'lazy': False,
         'minimal': False,
         'port_hide_changed': set([]),
         'posx': 246.25,
         'posy': -473.75,
         'priority': 0,
         'user_application': None},
   '__in__': {  'caption': 'In',
                'hide': True,
                'lazy': True,
                'minimal': False,
                'port_hide_changed': set([]),
                'posx': 140.0,
                'posy': -506.25,
                'priority': 0},
   '__out__': {  'caption': 'Out',
                 'hide': True,
                 'lazy': True,
                 'minimal': False,
                 'port_hide_changed': set([]),
                 'posx': 165.0,
                 'posy': -232.5,
                 'priority': 0}},
                             elt_value={7: [(1, 'None')], 8: [], 9: [], '__in__': [], '__out__': []},
                             lazy=True,
                             )




py_canview_plot_can = Factory(name='Plot Can File',
                description='Simple Plot of a can file',
                category='Visualisation.Caribu',
                nodemodule='py_canview',
                nodeclass='plot_can',
                inputs=[{'interface': IFileStr(filter="*.can", save=False), 'name': 'can file'}, {'interface': ISequence, 'name': 'colors', 'value': None}],
                outputs=[{'name': 'Canestra file'}],
                widgetmodule=None,
                widgetclass=None,
               )




