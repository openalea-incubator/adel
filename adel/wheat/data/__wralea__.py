# This file has been generated at Thu Feb 19 10:37:22 2009

from openalea.core import DataFactory as DF
from openalea.core import *

__name__ = 'alinea.adel.wheat.data'

__editable__ = True
__description__ = 'Wheat leaves measurement'
__license__ = ''
__url__ = ''
__alias__ = ['adel.wheat.data']
__version__ = '0.0.1'
__authors__ = 'Jessica Bertheloot'
__institutes__ = 'INRA'
__icon__ = ''

__all__ = ['_42296496', 'ToAdelParSo0N', '_42296528', '_42296656',
           'ToAdelParCa0N', '_42296720', '_42296240', '_114262320', '_42296592',
           '_42296688', '_42296624']

_42296496 = DF(uid="102d14e44f2011e6b469d4bed973e64a",
               name='ToAdelParSoN+',
               description='Input parameters for AdelWheat (cultivar Soissons under high N treatment)',
               editors=None,
               includes=None,
               )

ToAdelParSo0N = DF(uid="102d14e54f2011e6b469d4bed973e64a",
                   name='ToAdelParSo0N',
                   description='Input parameters for AdelWheat (cultivar Soissons under low N treatment)',
                   editors=None,
                   includes=None,
                   )

_42296528 = DF(uid="102d14e64f2011e6b469d4bed973e64a",
               name='ToAdelParCaN+',
               description='Input parameters for AdelWheat (cultivar Caphorn under high N treatment)',
               editors=None,
               includes=None,
               )

_42296656 = DF(uid="102d14e74f2011e6b469d4bed973e64a",
               name='AleaStr190.txt',
               description='One plant Adel Septo String',
               editors=None,
               includes=None,
               )

ToAdelParCa0N = DF(uid="1169117c4f2011e6b469d4bed973e64a",
                   name='ToAdelParCa0N',
                   description='Input parameters for AdelWheat (cultivar Caphorn under low N treatment) ',
                   editors=None,
                   includes=None,
                   )

_42296720 = DF(uid="1169117d4f2011e6b469d4bed973e64a",
               name='SampleString.str',
               description='Sample Lsytem string with LeafElement and StemElement',
               editors=None,
               includes=None,
               )

_42296240 = DF(uid="1169117e4f2011e6b469d4bed973e64a",
               name='testStr.txt',
               description='',
               editors=None,
               includes=None,
               )

_114262320 = DF(uid="1169117f4f2011e6b469d4bed973e64a",
                name='SRCa.RData',
                description='SR data for CapHorn (manip Jessica 2006)',
                editors=None,
                includes=None,
                )

_42296592 = DF(uid="116911804f2011e6b469d4bed973e64a",
               name='So99.RData',
               description='A set of X, Y coordinate of leaf midrib',
               editors=None,
               includes=None,
               )

_42296688 = DF(uid="116911814f2011e6b469d4bed973e64a",
               name='SRSo.RData',
               description='Curvilinear abscisse, radius of leaves ',
               editors=None,
               includes=None,
               )

_42296624 = DF(uid="116911824f2011e6b469d4bed973e64a",
               name='leaves_wheat.db',
               description='',
               editors=None,
               includes=None,
               )
