
# This file has been generated at Thu Apr 30 14:12:27 2009

from openalea.core import *


__name__ = 'alinea.adel.simulation'

__editable__ = True
__description__ = 'Simulation of a plant based on measurement using LSystem'
__license__ = 'CECILL'
__url__ = ''
__alias__ = ['adel.simulation']
__version__ = '0.0.1'
__authors__ = 'C. Fournier, C. Pradal'
__institutes__ = 'INRA, CIRAD, INRIA'
__icon__ = ''


__all__ = ['simulation_simulation', 'simulation_RunAdel','simulation_genString']



simulation_simulation = Factory(name='simulation',
                description='cpfg ',
                category='scene',
                nodemodule='simulation',
                nodeclass='simulation',
                inputs=[{'interface': IFileStr, 'name': 'lsystem', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'nb_plants', 'value': 1, 'desc': ''}],
                outputs=[{'interface': ITextStr, 'name': 'axial_tree', 'desc': 'Simulated plant'}],
                widgetmodule=None,
                widgetclass=None,
               )




simulation_RunAdel = Factory(name='RunAdel',
                description='',
                category='simulation',
                nodemodule='simulation',
                nodeclass='RunAdel',
                inputs=[{'interface': ISequence, 'name': 'dates', 'value': 1000, 'desc': ''},
                        {'interface': None, 'name': 'Parameter Rlist', 'desc': ''},
                        {'interface': IDict, 'name': 'Global parameters', 'value': {'senescence_leaf_shrink' : 0.5,'startLeaf' : -0.4, 'endLeaf' : 1.6, 'stemLeaf' : 1.2,'epsillon' : 1e-6}, 'desc': ''}
                        ],
                outputs=[{'interface': None, 'name': 'Lstrings', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )

simulation_genString = Factory(name='genString',
                description='',
                category='simulation',
                nodemodule='simulation',
                nodeclass='genString',
                inputs=[{'interface': None, 'name': 'R Canopy Table', 'desc': 'R Dataframe representing the canopy'}],
                outputs=[{'interface': None, 'name': 'Lstrings', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )



