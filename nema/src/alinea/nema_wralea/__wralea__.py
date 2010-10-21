
# This file has been generated at Fri Aug 13 14:50:32 2010

from openalea.core import *


__name__ = 'alinea.nema'

__editable__ = True
__description__ = 'Port to python/visualea of the nitrogen distribution model Nema'
__license__ = None
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = None
__authors__ = ''
__institutes__ = 'INRA, INRIA'
__icon__ = 'nema2520.png'


__all__ = ['nemafile_nemafile', 'plot_plotOut']



nemafile_nemafile = Factory(name='nemafile',
                description='run nema with files inputs',
                category='model',
                nodemodule='nemafile',
                nodeclass='nemafile',
                inputs=[{'interface': IDirStr, 'name': 'parameters', 'value': None, 'desc': ''}, {'interface': IDirStr, 'name': 'DrivingVariables', 'value': None, 'desc': ''}, {'interface': IDirStr, 'name': 'State0dd', 'value': None, 'desc': ''}, {'interface': IDirStr, 'name': 'outdir', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'Steps', 'value': 49, 'desc': 'Number of steps for the simulation'}, {'interface': IInt, 'name': 'substeps', 'value': 4, 'desc': 'number of substeps per step for the integration'}],
                outputs=[{'interface': IDirStr, 'name': 'outdir', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




plot_plotOut = Factory(name='plotOut',
                description='launch R function that plot Nema results',
                category='model',
                nodemodule='plot',
                nodeclass='plotOut',
                inputs=[{'interface': IDirStr, 'name': 'outdir', 'value': None, 'desc': 'directory that will contains the results'}],
                outputs=[{'interface': IDirStr, 'name': 'outdir', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




