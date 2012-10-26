
# This file has been generated at Fri Oct 26 18:01:13 2012

from openalea.core import *


__name__ = 'alinea.adel.plantgen'

__editable__ = True
__description__ = ''
__license__ = ''
__url__ = ''
__alias__ = ['plantgen']
__version__ = ''
__authors__ = ''
__institutes__ = ''
__icon__ = ''


__all__ = ['generate_adel_input_data_generate_adel_input_data']


generate_adel_input_data_generate_adel_input_data = Factory(name='generate_adel_input_data',
                authors='C. Chambon, M. Abichou and B. Andrieu',
                description='Generate ADEL input data.',
                category='data processing',
                nodemodule='generate_adel_input_data',
                nodeclass='generate_adel_input_data',
                inputs=None,
                outputs=({'interface': IInt, 'name': 'axis_table', 'desc': 'the axis table'}, {'interface': IInt, 'name': 'absolute_phen_table', 'desc': 'the absolute phen table'}, {'interface': IInt, 'name': 'relative_phen_table', 'desc': 'the relative phen table'}, {'interface': IInt, 'name': 'absolute_organ_dimensions_table', 'desc': 'the absolute dim table'}, {'interface': IInt, 'name': 'leaf_dynamic_parameters_table', 'desc': 'the leaf_dynamic_parameters table'}, {'interface': IInt, 'name': 'first_leaf_phen_table', 'desc': 'the first leaf phen table'}, {'interface': IInt, 'name': 'HS_GL_SSI_dynamic_table', 'desc': 'the HS_GL_SSI dynamic table'}, {'interface': IInt, 'name': 'relative_organ_dimensions_table', 'desc': 'the relative dim table'}, {'interface': IInt, 'name': 'tillering_dynamic_table', 'desc': 'the tillering dynamic table'}),
                widgetmodule=None,
                widgetclass=None,
               )



