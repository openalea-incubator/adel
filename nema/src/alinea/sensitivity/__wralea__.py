
# This file has been generated at Thu Aug 19 15:57:25 2010

from openalea.core import *


__name__ = 'alinea.nema.sensitivity'

__editable__ = True
__description__ = ''
__license__ = ''
__url__ = ''
__alias__ = []
__version__ = ''
__authors__ = ''
__institutes__ = ''
__icon__ = ''


__all__ = ['SensitivityTest_SensitivityTest', 'SensitivityAnalysisPB_SensitivityAnalysisPB']



SensitivityTest_SensitivityTest = Factory(name='SensitivityTest',
                description='',
                category='model',
                nodemodule='SensitivityTest',
                nodeclass='SensitivityTest',
                inputs=[{'interface': IDirStr, 'name': 'outdir', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDirStr, 'name': 'outdir', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




SensitivityAnalysisPB_SensitivityAnalysisPB = Factory(name='SensitivityAnalysisPB',
                description='Plackett?Burman technique for sensitivity analysis - initial test',
                category='model',
                nodemodule='SensitivityAnalysisPB',
                nodeclass='SensitivityAnalysisPB',
                inputs=[{'interface': IFileStr, 'name': 'inputfile', 'value': None, 'desc': ''}, {'interface': IDirStr, 'name': 'outdir', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDirStr, 'name': 'out', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




