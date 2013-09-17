
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


__all__ = ['plantgen_MIN_plantgen_MIN', 'plantgen_SHORT_plantgen_SHORT', 'plantgen_FULL_plantgen_FULL', 'read_plantgen_inputs_read_plantgen_inputs']


plantgen_MIN_plantgen_MIN = Factory(name='plantgen_MIN',
                authors='C. Chambon, M. Abichou and B. Andrieu',
                category='data processing',
                nodemodule='plantgen',
                nodeclass='gen_adel_input_data_from_min',
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

plantgen_SHORT_plantgen_SHORT = Factory(name='plantgen_SHORT',
                authors='C. Chambon, M. Abichou and B. Andrieu',
                category='data processing',
                nodemodule='plantgen',
                nodeclass='gen_adel_input_data_from_short',
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

plantgen_FULL_plantgen_FULL = Factory(name='plantgen_FULL',
                authors='C. Chambon, M. Abichou and B. Andrieu',
                category='data processing',
                nodemodule='plantgen',
                nodeclass='gen_adel_input_data_from_full',
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
                inputs=({'interface': IFileStr, 'name': 'inputs_filepath', 'desc': 'the Python module which describes the inputs of plantgen'}, 
                        {'interface': IStr, 'name': 'dynT_user_completeness', 'desc': 'the completeness of the table dynT_user'}),
                outputs=({'interface': IInterface, 'name': 'dynT_user'},
                         {'interface': IInterface, 'name': 'dimT_user'},
                         {'interface': IInt, 'name': 'plants_number'},
                         {'interface': IInt, 'name': 'plants_density'},
                         {'interface': IDict, 'name': 'decide_child_axis_probabilities'},
                         {'interface': IDict, 'name': 'MS_leaves_number_probabilities'},
                         {'interface': IInt, 'name': 'ears_density'},
                         {'interface': IDict, 'name': 'GL_number'},
                         {'interface': IFloat, 'name': 'delais_TT_stop_del_axis'},
                         {'interface': IFloat, 'name': 'TT_col_break'}),
                widgetmodule=None,
                widgetclass=None,
               )

