
# This file has been generated at Tue Feb 12 14:41:01 2008

from openalea.core import *


__name__ = 'alinea.adel.povray'
__alias__ = []

__editable__ = True
__version__ = '0.9.1'
__description__ = 'Use of povray to compute coverage rate'
__license__ = 'CECILL'
__authors__ = 'C. Pradal'
__url__ = ''
__institutes__ = 'CIRAD, INRIA'
 

__all__ = []



fit_fit = Factory(name='fit', 
                description='fit a leaf (x, y, s, r)', 
                category='fitting', 
                nodemodule='fit',
                nodeclass='fit',
                inputs=[{'interface': ISequence, 'name': 'x', 'value': [], 'desc': 'midrib axis'}, {'interface': ISequence, 'name': 'y', 'value': [], 'desc': 'midrib ordinate'}, {'interface': ISequence, 'name': 's', 'value': [], 'desc': 'midrib curvilinear abscisse'}, {'interface': ISequence, 'name': 'r', 'value': [], 'desc': 'midrib radius'}, {'interface': IInt, 'name': 'nb_points', 'value': 3, 'desc': 'target number of points after simplification'}],
                outputs=[{'interface': ISequence, 'name': 'leaf'}],
                )




fit_fit_leaf = Factory(name='fit midrib', 
                description='Fit a midrib with a spline curve of degree 5.', 
                category='fitting', 
                nodemodule='fit',
                nodeclass='fit_leaf',
                inputs=[{'interface': ISequence, 'name': 'leaf'}],
                outputs=[{'interface': ISequence, 'name': 'leaf'}, {'interface': IFloat, 'name': 'surface'}],
                widgetmodule=None,
                widgetclass=None,
                )




fit_simplify = Factory(name='simplification', 
                description='leaf simplification', 
                category='fitting', 
                nodemodule='fit',
                nodeclass='simplify',
                inputs=[{'interface': ISequence, 'name': 'leaf', 'value': [], 'desc': 'midrib (x, y, s, r)'}, {'interface': IInt, 'name': 'nb_points', 'value': 3, 'desc': 'target number of points after simplification'}],
                outputs=[{'interface': ISequence, 'name': 'leaf'}],
                widgetmodule=None,
                widgetclass=None,
                )




fit_leaves = Factory(name='fit leaves', 
                description='Replace measured midrib by fit and simplified curves', 
                category='fitting', 
                nodemodule='fit',
                nodeclass='fit_leaves',
                inputs=[{'interface': IDict, 'name': 'leaves', 'value': {}, 'desc': ''}, {'interface': IInt, 'name': 'nb_points', 'value': 7, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'leaves', 'desc': 'New database'}],
                widgetmodule=None,
                widgetclass=None,
                )




