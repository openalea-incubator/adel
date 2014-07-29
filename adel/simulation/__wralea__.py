
# This file has been generated at Thu Feb 17 22:34:17 2011

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


__all__ = ['simulation_simulation', 'simulation_RunAdel', 'AdelRunOptions_AdelRunOptions', 'simulation_genString']



simulation_simulation = Factory(name='simulation',
                authors='C. Fournier, C. Pradal (wralea authors)',
                description='cpfg ',
                category='scene',
                nodemodule='alinea.adel.simulation.simulation',
                nodeclass='simulation',
                inputs=[{'interface': IFileStr, 'name': 'lsystem', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'nb_plants', 'value': 1, 'desc': ''}],
                outputs=[{'interface': ITextStr, 'name': 'axial_tree', 'desc': 'Simulated plant'}],
                widgetmodule=None,
                widgetclass=None,
               )




simulation_RunAdel = Factory(name='RunAdel',
                authors='C. Fournier, C. Pradal (wralea authors)',
                description='',
                category='simulation',
                nodemodule='alinea.adel.simulation.simulation',
                nodeclass='RunAdel',
                inputs=[{'interface': ISequence, 'name': 'dates', 'value': 1000, 'desc': ''}, {'interface': None, 'name': 'Parameter Rlist', 'desc': ''}, {'interface': IDict, 'name': 'Global parameters', 'value': {'endLeaf': 1.6000000000000001,'endLeaf1': 1.6000000000000001, 'stemLeaf': 1.2, 'epsillon': 9.9999999999999995e-07, 'startLeaf': -0.40000000000000002, 'senescence_leaf_shrink': 0.5, 'HSstart_inclination_tiller': 1, 'rate_inclination_tiller': 30}, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'Lstrings', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




AdelRunOptions_AdelRunOptions = Factory(name='AdelRunOptions',
                authors='C. Fournier, C. Pradal (wralea authors)',
                description='Node interface for Adel General Option dict',
                category='parameter definition',
                nodemodule='alinea.adel.simulation.AdelRunOptions',
                nodeclass='AdelRunOptions',
                inputs=[{'interface': IFloat, 'name': 'startLeaf', 'value': -0.40000000000000002, 'desc': ''}, 
                        {'interface': IFloat, 'name': 'endLeaf', 'value': 1.6000000000000001, 'desc': ''}, 
                        {'interface': IFloat, 'name': 'stemLeaf', 'value': 1.2, 'desc': ''},
                        {'interface': IFloat, 'name': 'senescence_leaf_shrink', 'value': 0.5, 'desc': ''}, 
                        {'interface': IFloat, 'name': 'epsillon', 'value': (0, 1), 'desc': ''}, 
                        {'interface': IFloat, 'name': 'HSstart_inclination_tiller', 'value': 1, 'desc': ''}, 
                        {'interface': IFloat, 'name': 'rate_inclination_tiller', 'value': 30, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'Adel general option dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




simulation_genString = Factory(name='genString',
                authors='C. Fournier, C. Pradal (wralea authors)',
                description='',
                category='simulation',
                nodemodule='alinea.adel.simulation.simulation',
                nodeclass='genString',
                inputs=[{'interface': None, 'name': 'R Canopy Table', 'desc': 'R Dataframe representing the canopy'}],
                outputs=[{'interface': None, 'name': 'Lstrings', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




