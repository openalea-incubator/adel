
# This file has been generated at Wed Aug 25 20:11:15 2010

from openalea.core import *


__name__ = 'alinea.adel.wheat'

__editable__ = True
__description__ = '3D Reconstruction of a wheat canopy'
__license__ = ''
__url__ = ''
__alias__ = ['adel.wheat']
__version__ = '0.0.1'
__authors__ = 'Jessica Bertheloot'
__institutes__ = 'INRA'
__icon__ = ''


__all__ = ['Write_Table_Write_Table', 'extract_wheat_extract_leaf_info', 'selectOutput_Aggregate_selectOutput_Aggregate', 'CanMTGInterpreter_CanMTGInterpreter', 'Sum_EnergyDensities_Sum_EnergyDensities', 'Group_and_Apply_Group_and_Apply', 'To_Aggregation_Table_To_Aggregation_Table']



Write_Table_Write_Table = Factory(name='Write_Table',
                description='Write a table (dict) in a file',
                category='io',
                nodemodule='Write_Table',
                nodeclass='Write_Table',
                inputs=[{'interface': IData, 'name': 'Table', 'value': None, 'desc': ''},
                        {'interface': IFileStr, 'name': 'File name', 'value': None, 'desc': ''},
                        {'interface': IBool, 'name': 'first time called ?', 'value': True, 'desc': 'determine header writing and mode for opening file (w,a)'}],
                outputs=[{'interface': IFileStr, 'name': 'txtFile', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




extract_wheat_extract_leaf_info = Factory(name='extract_leaves',
                description='Build a dict of leaves (x, y, s, r)  from R data files',
                category='Unclassified',
                nodemodule='extract_wheat',
                nodeclass='extract_leaf_info',
                inputs=[{'interface': IFileStr, 'name': 'xy', 'value': None, 'desc': ''}, {'interface': IFileStr, 'name': 'sr', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'leaves', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




selectOutput_Aggregate_selectOutput_Aggregate = Factory(name='selectOutput_Aggregate',
                description='Selection of one output of the aggregate dictionnary according to the keys and conversion into a sequence',
                category='io',
                nodemodule='selectOutput_Aggregate',
                nodeclass='selectOutput_Aggregate',
                inputs=[{'interface': IDict, 'name': 'table', 'value': None, 'desc': ''}, {'interface': IEnumStr(enum=['Plant', 'LeafElement', 'Metamer', 'StemElement', 'Axe']), 'name': 'variable', 'value': None, 'desc': ''}],
                outputs=[{'interface': ISequence, 'name': 'Variable Values', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




CanMTGInterpreter_CanMTGInterpreter = Factory(name='CanMTGInterpreter',
                description='Generates a CanMTG from unctions and Lsystem string',
                category='scene.PGL Object Generator',
                nodemodule='CanMTGInterpreter',
                nodeclass='CanMTGInterpreter',
                inputs=[{'interface': IData, 'name': 'Symbols', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'LSystem String', 'value': None, 'desc': ''}],
                outputs=[{'interface': IData, 'name': 'CanMTG', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




Sum_EnergyDensities_Sum_EnergyDensities = Factory(name='Sum_EnergyDensities',
                description='Return the mean energy density for two types of elements',
                category='data processing',
                nodemodule='Sum_EnergyDensities',
                nodeclass='Sum_EnergyDensities',
                inputs=[{'interface': ISequence, 'name': 'energy1', 'value': None, 'desc': ''}, {'interface': ISequence, 'name': 'area1', 'value': None, 'desc': ''}, {'interface': ISequence, 'name': 'energy2', 'value': None, 'desc': ''}, {'interface': ISequence, 'name': 'area2', 'value': None, 'desc': ''}],
                outputs=[{'interface': IFloat, 'name': 'Energy density', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




Group_and_Apply_Group_and_Apply = Factory(name='Group and Apply',
                description='',
                category='Unclassified',
                nodemodule='Group_and_Apply',
                nodeclass='Group_and_Apply',
                inputs=[{'interface': IDict, 'name': 'table', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'properties', 'value': None, 'desc': ''}, {'interface': IFunction, 'name': 'function', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'data', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'OUT1', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




To_Aggregation_Table_To_Aggregation_Table = Factory(name='To_Aggregation_Table',
                description='Return a dictionnary from the canMTG',
                category='io',
                nodemodule='To_Aggregation_Table',
                nodeclass='To_Aggregation_Table',
                inputs=[{'interface': IData, 'name': 'CanMTG', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'Table', 'desc': ''},
                         {'name': 'CanMTG', 'desc': 'MTG with new properties'}],
                widgetmodule=None,
                widgetclass=None,
               )




