
# This file has been generated at Fri Oct 26 18:01:13 2012

from openalea.core import *


__name__ = 'alinea.adel.plantgen'

__editable__ = True
__description__ = ''
__license__ = ''
__url__ = ''
__alias__ = ['plantgen']
__version__ = ''
__authors__ = 'M. Abichou, B. Andrieu, C. Chambon'
__institutes__ = ''
__icon__ = ''


__all__ = ['plantgen_plantgen', 'read_plantgen_inputs_read_plantgen_inputs', 'plantgen2adel_plantgen2adel']


plantgen_plantgen = Factory(name='plantgen',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                category='data processing',
                nodemodule='alinea.adel.plantgen.plantgen_interface',
                nodeclass='gen_adel_input_data',
                inputs=None,
                outputs=({'interface': None, 'name': 'axeT', 'desc': 'the axeT dataframe'}, 
                         {'interface': None, 'name': 'dimT', 'desc': 'the dimT dataframe'}, 
                         {'interface': None, 'name': 'phenT', 'desc': 'the phenT dataframe'}, 
                         {'interface': None, 'name': 'phenT_abs', 'desc': 'the phenT_abs dataframe'}, 
                         {'interface': None, 'name': 'dynT', 'desc': 'the dynT dataframe'}, 
                         {'interface': None, 'name': 'phenT_first', 'desc': 'the phenT_first dataframe'}, 
                         {'interface': None, 'name': 'HS_GL_SSI_T', 'desc': 'the HS_GL_SSI_T dataframe'}, 
                         {'interface': None, 'name': 'tilleringT', 'desc': 'the tilleringT dataframe'},
                         {'interface': None, 'name': 'cardinalityT', 'desc': 'the cardinalityT dataframe'},
                         {'interface': IDict, 'name': 'config', 'desc': 'the configuration used for the construction'}),
                widgetmodule=None,
                widgetclass=None,
               )

read_plantgen_inputs_read_plantgen_inputs = Factory(name='read_plantgen_inputs',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                category='data processing',
                nodemodule='alinea.adel.plantgen.plantgen_interface',
                nodeclass='read_plantgen_inputs',
                inputs=({'interface': IFileStr, 'name': 'inputs_filepath', 'desc': 'the Python module which describes the inputs of plantgen'},
                        {'interface': IFileStr, 'name': 'dynT_user_filepath', 'desc': 'the file path of the leaf dynamic parameters set by the user.'},
                        {'interface': IFileStr, 'name': 'dimT_user_filepath', 'desc': 'the file path of the dimensions of the organs set by the user.'}),
                outputs=({'interface': IInterface, 'name': 'dynT_user'},
                         {'interface': IInterface, 'name': 'dimT_user'},
                         {'interface': IInt, 'name': 'plants_number'},
                         {'interface': IInt, 'name': 'plants_density'},
                         {'interface': IDict, 'name': 'decide_child_axis_probabilities'},
                         {'interface': IDict, 'name': 'MS_leaves_number_probabilities'},
                         {'interface': IInt, 'name': 'ears_density'},
                         {'interface': IDict, 'name': 'GL_number'},
                         {'interface': IFloat, 'name': 'delais_TT_stop_del_axis'},
                         {'interface': IFloat, 'name': 'TT_col_break'},
                         {'interface': IDict, 'name': 'inner_params'}),
                widgetmodule=None,
                widgetclass=None,
               )

plantgen2adel_plantgen2adel = Factory(name='plantgen2adel',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                category='data processing',
                nodemodule='alinea.adel.plantgen.plantgen_interface',
                nodeclass='plantgen2adel',
                inputs=({'interface': None, 'name': 'axeT', 'desc': 'the axeT dataframe'},
                        {'interface': None, 'name': 'dimT', 'desc': 'the dimT dataframe'},
                        {'interface': None, 'name': 'phenT', 'desc': 'the phenT dataframe'}),
                outputs=({'interface': None, 'name': 'axeT_adel', 'desc': 'the axeT dataframe in adel-like format'},
                         {'interface': None, 'name': 'dimT_adel', 'desc': 'the dimT dataframe in adel-like format'},
                         {'interface': None, 'name': 'phenT_adel', 'desc': 'the phenT dataframe in adel-like format'}),
                widgetmodule=None,
                widgetclass=None,
               )

