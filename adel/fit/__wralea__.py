
# This file has been generated at Mon Nov 07 12:25:28 2011

from openalea.core import *


__name__ = 'alinea.adel.fitting'

__editable__ = True
__description__ = 'midrib fitting and simplification'
__license__ = 'CECILL'
__url__ = ''
__alias__ = ['adel.fitting']
__version__ = '0.0.1'
__authors__ = 'C. Pradal, C. Fournier'
__institutes__ = 'INRA, CIRAD, INRIA'
__icon__ = ''


__all__ = ['fit_fit_leaves', 'fit_fit_leaf', 'fit_fit', 'dimension_fitting_dimension_fitting', 'fit_simplify']



fit_fit_leaves = Factory(name='fit leaves',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Replace measured midrib by fit and simplified curves',
                category='fitting',
                nodemodule='fit',
                nodeclass='fit_leaves',
                inputs=[{'interface': IDict, 'name': 'leaves', 'value': {}, 'desc': ''}, {'interface': IInt, 'name': 'nb_points', 'value': 7, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'leaves', 'desc': 'New database'}],
                widgetmodule=None,
                widgetclass=None,
               )




fit_fit_leaf = Factory(name='fit midrib',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Fit a midrib with a spline curve of degree 5.',
                category='fitting',
                nodemodule='fit',
                nodeclass='fit_leaf',
                inputs=[{'interface': ISequence, 'name': 'leaf'}],
                outputs=[{'interface': ISequence, 'name': 'leaf'}, {'interface': IFloat, 'name': 'surface'}],
                widgetmodule=None,
                widgetclass=None,
               )




fit_fit = Factory(name='fit',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='fit a leaf (x, y, s, r)',
                category='fitting',
                nodemodule='fit',
                nodeclass='fit',
                inputs=[{'interface': ISequence, 'name': 'x', 'value': [], 'desc': 'midrib axis'}, {'interface': ISequence, 'name': 'y', 'value': [], 'desc': 'midrib ordinate'}, {'interface': ISequence, 'name': 's', 'value': [], 'desc': 'midrib curvilinear abscisse'}, {'interface': ISequence, 'name': 'r', 'value': [], 'desc': 'midrib radius'}, {'interface': IInt, 'name': 'nb_points', 'value': 3, 'desc': 'target number of points after simplification'}],
                outputs=[{'interface': ISequence, 'name': 'leaf'}],
                widgetmodule=None,
                widgetclass=None,
               )




dimension_fitting_dimension_fitting = Factory(name='dimension_fitting',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='dimension_fitting',
                nodeclass='dimension_fitting',
                inputs=[{'interface': None, 'name': 'ref_axis', 'value': None, 'desc': ''}, {'interface': None, 'name': 'dimT_tofit', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dimT_fitted_df', 'desc': 'pandas.dataframe'}],
                widgetmodule=None,
                widgetclass=None,
               )




fit_simplify = Factory(name='simplification',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='leaf simplification',
                category='fitting',
                nodemodule='fit',
                nodeclass='simplify',
                inputs=[{'interface': ISequence, 'name': 'leaf', 'value': [], 'desc': 'midrib (x, y, s, r)'}, {'interface': IInt, 'name': 'nb_points', 'value': 3, 'desc': 'target number of points after simplification'}],
                outputs=[{'interface': ISequence, 'name': 'leaf'}],
                widgetmodule=None,
                widgetclass=None,
               )




