# This file has been generated at Tue Dec 15 12:18:16 2009

from openalea.core import Factory as Fa
from openalea.core import ITextStr

__name__ = 'alinea.adel'

__editable__ = True
__description__ = 'Simulation of plant development and reconstruction'
__license__ = 'CECILL'
__url__ = 'http://openalea.gforge.inria.fr/'
__alias__ = ['adel']
__version__ = '0.0.1'
__authors__ = 'C. Fournier, C. Pradal'
__institutes__ = 'INRA, CIRAD'
__icon__ = 'adel.png'

__all__ = ['lpy2mtg_lpy2mtg']

lpy2mtg_lpy2mtg = Fa(uid="45880b624f1f11e6b469d4bed973e64a",
                     name='lpy2mtg',
                     authors='C. Fournier and C. Pradal (wralea authors)',
                     description='',
                     category='Unclassified',
                     nodemodule='alinea.adel.lpy2mtg',
                     nodeclass='lpy2mtg',
                     inputs=[{'interface': ITextStr, 'name': 'AxialTree',
                              'value': None, 'desc': ''},
                             {'interface': None, 'name': 'LSystem',
                              'value': None, 'desc': ''}],
                     outputs=[{'interface': None, 'name': 'mtg', 'desc': ''}],
                     widgetmodule=None,
                     widgetclass=None,
                     )
