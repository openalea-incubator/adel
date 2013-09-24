
# This file has been generated at Fri Oct 26 18:01:13 2012

from openalea.core import *


__name__ = 'alinea.adel.plantgen'

__editable__ = True
__description__ = ''
__license__ = ''
__url__ = ''
__alias__ = ['plantgen']
__version__ = ''
__authors__ = 'C. Chambon, M. Abichou and B. Andrieu'
__institutes__ = ''
__icon__ = ''


__all__ = ['plantgen_plantgen', 'read_plantgen_inputs_read_plantgen_inputs']


plantgen_plantgen = Factory(name='plantgen',
                authors='C. Chambon, M. Abichou and B. Andrieu',
                category='data processing',
                nodemodule='plantgen',
                nodeclass='gen_adel_input_data',
                inputs=None,
                outputs=({'interface': None, 'name': 'axeT', 'desc': 'the axeT dataframe'}, 
                         {'interface': None, 'name': 'dimT', 'desc': 'the dimT dataframe'}, 
                         {'interface': None, 'name': 'phenT', 'desc': 'the phenT dataframe'}, 
                         {'interface': None, 'name': 'phenT_abs', 'desc': 'the phenT_abs dataframe'}, 
                         {'interface': None, 'name': 'dimT_abs', 'desc': 'the dimT_abs dataframe'}, 
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
                authors='C. Chambon, M. Abichou and B. Andrieu',
                category='data processing',
                nodemodule='plantgen',
                nodeclass='read_plantgen_inputs',
                inputs=({'interface': IFileStr, 'name': 'inputs_filepath', 'desc': 'the Python module which describes the inputs of plantgen'},),
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
                         {'interface': IStr, 'name': 'dynT_user_completeness'},
                         {'interface': IStr, 'name': 'dimT_user_completeness'},
                         {'interface': IDict, 'name': 'inner_params'}),
                widgetmodule=None,
                widgetclass=None,
               )

