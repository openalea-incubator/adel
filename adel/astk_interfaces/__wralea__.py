
# This file has been generated at Wed Nov 07 14:27:04 2012

from openalea.core import *


__name__ = 'alinea.adel.astk_interfaces'

__editable__ = True
__description__ = 'Interfaces complying to astk plant_interface'
__license__ = ''
__url__ = ''
__version__ = '0.0.1'
__authors__ = 'Christian Fournier'
__institutes__ = 'INRA'
__icon__ = ''


__all__ = []


astk_AdelWheat = Factory(name='AdelWheat',
                nodemodule='alinea.adel.astk_interface',
                nodeclass='adelwheat_node',
               )
__all__.append('astk_AdelWheat')

