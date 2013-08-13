
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


__all__ = ['plantgen_MIN_plantgen_MIN', 'plantgen_SHORT_plantgen_SHORT', 'plantgen_FULL_plantgen_FULL']


plantgen_MIN_plantgen_MIN = Factory(name='plantgen_MIN',
                authors='C. Chambon, M. Abichou and B. Andrieu',
                category='data processing',
                nodemodule='plantgen',
                nodeclass='gen_adel_input_data_from_min',
                inputs=None,
                outputs=({'interface': IInt, 'name': 'axeT', 'desc': 'the axeT dataframe'}, 
                         {'interface': IInt, 'name': 'dimT', 'desc': 'the dimT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT', 'desc': 'the phenT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT_abs', 'desc': 'the phenT_abs dataframe'}, 
                         {'interface': IInt, 'name': 'dimT_abs', 'desc': 'the dimT_abs dataframe'}, 
                         {'interface': IInt, 'name': 'dynT', 'desc': 'the dynT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT_first', 'desc': 'the phenT_first dataframe'}, 
                         {'interface': IInt, 'name': 'HS_GL_SSI_T', 'desc': 'the HS_GL_SSI_T dataframe'}, 
                         {'interface': IInt, 'name': 'tilleringT', 'desc': 'the tilleringT dataframe'},
                         {'interface': IInt, 'name': 'cardinalityT', 'desc': 'the cardinalityT dataframe'}),
                widgetmodule=None,
                widgetclass=None,
               )

plantgen_SHORT_plantgen_SHORT = Factory(name='plantgen_SHORT',
                authors='C. Chambon, M. Abichou and B. Andrieu',
                category='data processing',
                nodemodule='plantgen',
                nodeclass='gen_adel_input_data_from_short',
                inputs=None,
                outputs=({'interface': IInt, 'name': 'axeT', 'desc': 'the axeT dataframe'}, 
                         {'interface': IInt, 'name': 'dimT', 'desc': 'the dimT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT', 'desc': 'the phenT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT_abs', 'desc': 'the phenT_abs dataframe'}, 
                         {'interface': IInt, 'name': 'dimT_abs', 'desc': 'the dimT_abs dataframe'}, 
                         {'interface': IInt, 'name': 'dynT', 'desc': 'the dynT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT_first', 'desc': 'the phenT_first dataframe'}, 
                         {'interface': IInt, 'name': 'HS_GL_SSI_T', 'desc': 'the HS_GL_SSI_T dataframe'}, 
                         {'interface': IInt, 'name': 'tilleringT', 'desc': 'the tilleringT dataframe'},
                         {'interface': IInt, 'name': 'cardinalityT', 'desc': 'the cardinalityT dataframe'}),
                widgetmodule=None,
                widgetclass=None,
               )

plantgen_FULL_plantgen_FULL = Factory(name='plantgen_FULL',
                authors='C. Chambon, M. Abichou and B. Andrieu',
                category='data processing',
                nodemodule='plantgen',
                nodeclass='gen_adel_input_data_from_full',
                inputs=None,
                outputs=({'interface': IInt, 'name': 'axeT', 'desc': 'the axeT dataframe'}, 
                         {'interface': IInt, 'name': 'dimT', 'desc': 'the dimT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT', 'desc': 'the phenT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT_abs', 'desc': 'the phenT_abs dataframe'}, 
                         {'interface': IInt, 'name': 'dimT_abs', 'desc': 'the dimT_abs dataframe'}, 
                         {'interface': IInt, 'name': 'dynT', 'desc': 'the dynT dataframe'}, 
                         {'interface': IInt, 'name': 'phenT_first', 'desc': 'the phenT_first dataframe'}, 
                         {'interface': IInt, 'name': 'HS_GL_SSI_T', 'desc': 'the HS_GL_SSI_T dataframe'}, 
                         {'interface': IInt, 'name': 'tilleringT', 'desc': 'the tilleringT dataframe'},
                         {'interface': IInt, 'name': 'cardinalityT', 'desc': 'the cardinalityT dataframe'}),
                widgetmodule=None,
                widgetclass=None,
               )



