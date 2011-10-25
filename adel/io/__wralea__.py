
# This file has been generated at Tue Oct 26 13:41:12 2010

from openalea.core import *


__name__ = 'alinea.adel.io'

__editable__ = True
__description__ = 'Manage input and output for adel such as meseasurement.'
__license__ = 'CECILL'
__url__ = ''
__alias__ = ['adel.io']
__version__ = '0.0.1'
__authors__ = 'C. Pradal, C. Fournier'
__institutes__ = 'INRA, CIRAD, INRIA'
__icon__ = ''


__all__ = ['io_dataframe', 'io_dataframeAsdict','io_csvAsDict', 'io_saveRData', 'io_RlistAsDict', 'io_to_canestra', 'mylist_mylist', 'dataFrameAsDict_dataFrameAsDict', 'io_readRData', 'GetAdelString_GetAdelString', 'io_load_leaf_data', 'io_to_plantgl','io_canL2canS']

io_canL2canS = Factory(name='canL2canS',
                description='transform canopy Table Rdata to Surface canopy Table',
                category='io',
                nodemodule='io',
                nodeclass='canL2canS',
                inputs=[{'interface': None, 'name': 'canopy table (Rdataframe)', 'desc': ''},{'interface': None, 'name': 'sr database (RList)', 'desc': ''},{'interface': IFloat, 'name': 'senescence leaf shrink', 'desc': '','value': 1}],
                outputs=[{'interface': None, 'name': 'dataframe', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )


io_dataframe = Factory(name='Rdataframe',
                description='returns a dataframe (rpy2 object) from a dictionary containing named vectors of values',
                category='io',
                nodemodule='io',
                nodeclass='dataframe',
                inputs=[{'interface': IDict, 'name': 'Dict', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dataframe', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )

io_dataframeAsdict = Factory(name='dataframeAsdict',
                description='returns a dictionary containing named vectors of values from adataframe (rpy2 object)',
                category='io',
                nodemodule='io',
                nodeclass='dataframeAsdict',
                inputs=[{'interface': IDict, 'name': 'dataframe', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )



io_csvAsDict = Factory(name='CsvAsDict',
                description='returns a dictionary containing the columns of the csv file',
                category='io',
                nodemodule='io',
                nodeclass='csvAsDict',
                inputs=[{'interface': IFileStr, 'name': 'Csv File name', 'desc': ''}, {'interface': IInt, 'name': 'csv type', 'value': 1, 'desc': '1 for coma separated/dot as decimal, 2 for ; separated, coma decimal'}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_saveRData = Factory(name='saveRData',
                description='',
                category='io',
                nodemodule='io',
                nodeclass='saveRData',
                inputs=[{'name': 'RObject', 'desc': 'R object'}, {'interface': IStr, 'name': 'Name of the savec object', 'value': 'Robj', 'desc': ''}, {'interface': IFileStr, 'name': 'RData file', 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'name', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_RlistAsDict = Factory(name='RlistAsDict',
                description='returns a dictionary containing the Robjects of the Rlist',
                category='io',
                nodemodule='io',
                nodeclass='RlistAsDict',
                inputs=[{'interface': IFileStr, 'name': 'R list', 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_to_canestra = Factory(name='to_canestra',
                description='Convert Canestra Scene to a CanFile format stream.',
                category='io',
                nodemodule='io',
                nodeclass='to_canestra',
                inputs=[{'interface': IDict, 'name': 'scene', 'value': {}, 'desc': 'CanestraScene'}],
                outputs=[{'interface': ITextStr, 'name': 'string', 'desc': 'Canestra File'}],
                widgetmodule=None,
                widgetclass=None,
               )




mylist_mylist = Factory(name='mylist',
                description='convert numpy array inot list',
                category='Unclassified',
                nodemodule='mylist',
                nodeclass='mylist',
                inputs=[{'interface': None, 'name': 'array', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'list', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




dataFrameAsDict_dataFrameAsDict = Factory(name='dataFrameAsDict',
                description='',
                category='Unclassified',
                nodemodule='dataFrameAsDict',
                nodeclass='dataFrameAsDict',
                inputs=[{'interface': None, 'name': 'RDataframe', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_readRData = Factory(name='readRData',
                description='return a dictionary containing the Robject of the RData file',
                category='io',
                nodemodule='io',
                nodeclass='readRData',
                inputs=[{'interface': IFileStr, 'name': 'RData file', 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




GetAdelString_GetAdelString = Factory(name='GetAdelString',
                description='Reurn the content of a numbered lsystem string file',
                category='io',
                nodemodule='GetAdelString',
                nodeclass='GetAdelString',
                inputs=[{'interface': IDirStr, 'name': 'LsysDir', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'IterNumber', 'value': None, 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'LsysString', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_load_leaf_data = Factory(name='load leaf data',
                description='Load leaf data obtained by measurement',
                category='io',
                nodemodule='io',
                nodeclass='load_leaf_data',
                inputs=[{'interface': IFileStr, 'name': 'filename', 'value': None, 'desc': 'filename of a pickle file'}],
                outputs=[{'interface': IDict, 'name': 'leaves', 'desc': 'A dictionary associating leaf rank to a list of leaf values (x, y, s, r)'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_to_plantgl = Factory(name='to_plantgl',
                description='Adapt a Canestra scene to a PlantGL scene',
                category='io',
                nodemodule='io',
                nodeclass='to_plantgl',
                inputs=[{'name': 'scene', 'desc': 'Canestra Scene'}, {'interface': IRGBColor, 'name': 'leaf_color', 'value': (0, 180, 0)}, {'interface': IRGBColor, 'name': 'stem_color', 'value': (0, 130, 0)}, {'interface': IRGBColor, 'name': 'soil_color', 'value': (170, 85, 0)}],
                outputs=[{'interface': IInterface, 'name': 'scene', 'desc': 'PlantGL scene'}],
                widgetmodule=None,
                widgetclass=None,
               )



params2mtg = Factory(name='mtg (params)',
                category='simulation',
                nodemodule='io',
                nodeclass='mtg_factory',
               )
__all__.append('params2mtg')

tht = Factory(name='update thermal time',
                category='simulation',
                nodemodule='io',
                nodeclass='thermal_time',
               )
__all__.append('tht')

duplicate = Factory(name='duplicate mtg',
                category='simulation',
                nodemodule='io',
                nodeclass='duplicate',
                inputs=[dict(name='mtg'), dict(name='nb_plants', interface='IInt', value=1)],
                outputs=[{'name': 'mtg', 'desc': 'duplicated mtg'}],
               )
__all__.append('duplicate')


lpy2mtg_lpy2mtg = Factory(name='lpy2mtg',
                description='aggregate lpy outputs into an mtg',
                category='data i/o',
                nodemodule='io',
                nodeclass='lpy2mtg',
                inputs=[{'interface': None, 'name': 'axial tree', 'value': None, 'desc': ''}, {'interface': None, 'name': 'lsystem', 'value': None, 'desc': ''}, {'interface': None, 'name': 'scene', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'mtg', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )


__all__.append('lpy2mtg_lpy2mtg')

mtg2axial_mtg2lpy = Factory(name='mtg2axial',
                description='Convert MTg to axial tree along with spec in lysystem',
                category='data i/o',
                nodemodule='io',
                nodeclass='mtg2lpy',
                inputs=[{'interface': None, 'name': 'mtg', 'value': None, 'desc': ''}, {'interface': None, 'name': 'lsystem', 'value': None, 'desc': ''}, {'interface': None, 'name': 'axial tree', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'axial tree', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )

__all__.append('mtg2axial_mtg2lpy')

