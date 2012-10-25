
# This file has been generated at Wed Apr 25 08:10:00 2012

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


__all__ = ['io_duplicate', 'io_csv2pandasDataframe', 'io_canL2canS', 'io_mtg_factory', 'io_csvAsDict', 'io_RlistAsDict', 'io_dataframeAsdict', 'io_dataframe', 'io_pandasDataframe2csv', 'io_thermal_time', 'io_to_canestra', 'io_saveRData', 'io_lpy2mtg', 'PairAsDict_PairAsDict', 'io_readRData', 'io_apply_property', 'GetAdelString_GetAdelString', 'io_load_leaf_data', 'io_mtg2lpy', 'io_to_plantgl']



io_duplicate = Factory(name='duplicate mtg',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='simulation',
                nodemodule='io',
                nodeclass='duplicate',
                inputs=[{'name': 'mtg'}, {'interface': 'IInt', 'name': 'nb_plants', 'value': 1}],
                outputs=[{'name': 'mtg', 'desc': 'duplicated mtg'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_csv2pandasDataframe = Factory(name='csv2pandasDataframe',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Read CSV (comma-separated) file into DataFrame.',
                category='data i/o',
                nodemodule='io',
                nodeclass='csv2pandasDataframe',
                inputs=[{'interface': IFileStr, 'name': 'csv_filepath', 'value': None, 'desc': 'The filepath of the csv file to import.'}, {'interface': ISequence, 'name': 'index_col', 'value': None, 'desc': 'Column to use as the row labels of the DataFrame. If a sequence is given, a MultiIndex is used.'}, {'interface': ISequence, 'name': 'na_values', 'value': None, 'desc': 'List of additional strings to recognize as NA/NaN.'}, {'interface': IBool, 'name': 'parse_dates', 'value': False, 'desc': 'Attempt to parse dates in the index column(s).'}],
                outputs=[{'interface': None, 'name': 'dataframe', 'desc': 'A pandas.DataFrame instance which represents the csv file.'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_canL2canS = Factory(name='canL2canS',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='transform canopy Table Rdata to Surface canopy Table',
                category='io',
                nodemodule='io',
                nodeclass='canL2canS',
                inputs=[{'interface': None, 'name': 'canopy table (Rdataframe)', 'desc': ''}, {'interface': None, 'name': 'sr database (RList)', 'desc': ''}, {'interface': IFloat, 'name': 'senescence leaf shrink', 'value': 1, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dataframe', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_mtg_factory = Factory(name='mtg (params)',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='simulation',
                nodemodule='io',
                nodeclass='mtg_factory',
                inputs=[{'interface': IDict, 'name': 'Canopy table'}, {'interface': IInt, 'name': 'number of sectors', 'value': 1}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




io_csvAsDict = Factory(name='CsvAsDict',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dictionary containing the columns of the csv file',
                category='io',
                nodemodule='io',
                nodeclass='csvAsDict',
                inputs=[{'interface': IFileStr, 'name': 'Csv File name', 'desc': ''}, {'interface': IInt, 'name': 'csv type', 'value': 1, 'desc': '1 for coma separated/dot as decimal, 2 for ; separated, coma decimal'}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_RlistAsDict = Factory(name='RlistAsDict',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dictionary containing the Robjects of the Rlist',
                category='io',
                nodemodule='io',
                nodeclass='RlistAsDict',
                inputs=[{'interface': IFileStr, 'name': 'R list', 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_dataframeAsdict = Factory(name='dataframeAsdict',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dictionary containing named vectors of values from adataframe (rpy2 object)',
                category='io',
                nodemodule='io',
                nodeclass='dataframeAsdict',
                inputs=[{'interface': IDict, 'name': 'dataframe', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_dataframe = Factory(name='Rdataframe',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dataframe (rpy2 object) from a dictionary containing named vectors of values',
                category='io',
                nodemodule='io',
                nodeclass='dataframe',
                inputs=[{'interface': IDict, 'name': 'Dict', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dataframe', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_pandasDataframe2csv = Factory(name='pandasDataframe2csv',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Write a DataFrame to a comma-separated values (csv) file.',
                category='data i/o',
                nodemodule='io',
                nodeclass='pandasDataframe2csv',
                inputs=[{'interface': None, 'name': 'dataframe', 'value': None, 'desc': 'The DataFrame to write.'}, {'interface': IFileStr, 'name': 'csv_filepath', 'value': None, 'desc': 'The file path where the Dataframe is written.'}, {'interface': IStr, 'name': 'na_rep', 'value': None, 'desc': 'Missing data replacement.'}, {'interface': IBool, 'name': 'index', 'value': True, 'desc': 'Write row names (index)'}, {'interface': ISequence, 'name': 'index_label', 'value': None, 'desc': 'Column label for index column(s) if desired. If None is given, and header and index are True, then the index names are used. A sequence should be given if the DataFrame uses MultiIndex.'}],
                outputs=[{'interface': IFileStr, 'name': 'csv_filepath', 'desc': 'The file path where the Dataframe is written.'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_thermal_time = Factory(name='update thermal time',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='simulation',
                nodemodule='io',
                nodeclass='thermal_time',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




io_to_canestra = Factory(name='to_canestra',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Convert Canestra Scene to a CanFile format stream.',
                category='io',
                nodemodule='io',
                nodeclass='to_canestra',
                inputs=[{'interface': IDict, 'name': 'scene', 'value': {}, 'desc': 'CanestraScene'}],
                outputs=[{'interface': ITextStr, 'name': 'string', 'desc': 'Canestra File'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_saveRData = Factory(name='saveRData',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='io',
                nodemodule='io',
                nodeclass='saveRData',
                inputs=[{'name': 'RObject', 'desc': 'R object'}, {'interface': IStr, 'name': 'Name of the savec object', 'value': 'Robj', 'desc': ''}, {'interface': IFileStr, 'name': 'RData file', 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'name', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_lpy2mtg = Factory(name='lpy2mtg',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='aggregate lpy outputs into an mtg',
                category='data i/o',
                nodemodule='io',
                nodeclass='lpy2mtg',
                inputs=[{'interface': None, 'name': 'axial tree', 'value': None, 'desc': ''}, {'interface': None, 'name': 'lsystem', 'value': None, 'desc': ''}, {'interface': None, 'name': 'scene', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'mtg', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




PairAsDict_PairAsDict = Factory(name='PairAsDict',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dict from a key,value pair',
                category='data i/o',
                nodemodule='PairAsDict',
                nodeclass='PairAsDict',
                inputs=[{'interface': ITuple, 'name': 'Pair', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_readRData = Factory(name='readRData',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='return a dictionary containing the Robject of the RData file',
                category='io',
                nodemodule='io',
                nodeclass='readRData',
                inputs=[{'interface': IFileStr, 'name': 'RData file', 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_apply_property = Factory(name='apply_property',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='data i/o, MTG',
                nodemodule='io',
                nodeclass='apply_property',
                inputs=[{'interface': None, 'name': 'mtg', 'value': None, 'desc': 'MTG'}, {'interface': IStr, 'name': 'property name', 'value': '', 'desc': 'A MTG property'}, {'interface': IFunction, 'name': 'function', 'desc': 'function to apply on the property'}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': 'the output property'}],
                widgetmodule=None,
                widgetclass=None,
               )




GetAdelString_GetAdelString = Factory(name='GetAdelString',
                authors='C. Pradal, C. Fournier (wralea authors)',
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
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Load leaf data obtained by measurement',
                category='io',
                nodemodule='io',
                nodeclass='load_leaf_data',
                inputs=[{'interface': IFileStr, 'name': 'filename', 'value': None, 'desc': 'filename of a pickle file'}],
                outputs=[{'interface': IDict, 'name': 'leaves', 'desc': 'A dictionary associating leaf rank to a list of leaf values (x, y, s, r)'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_mtg2lpy = Factory(name='mtg2axial',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Convert MTg to axial tree along with spec in lysystem',
                category='data i/o',
                nodemodule='io',
                nodeclass='mtg2lpy',
                inputs=[{'interface': None, 'name': 'mtg', 'value': None, 'desc': ''}, {'interface': None, 'name': 'lsystem', 'value': None, 'desc': ''}, {'interface': None, 'name': 'axial tree', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'axial tree', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_to_plantgl = Factory(name='to_plantgl',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Adapt a Canestra scene to a PlantGL scene',
                category='io',
                nodemodule='io',
                nodeclass='to_plantgl',
                inputs=[{'name': 'scene', 'desc': 'Canestra Scene'}, {'interface': IRGBColor, 'name': 'leaf_color', 'value': (0, 180, 0)}, {'interface': IRGBColor, 'name': 'stem_color', 'value': (0, 130, 0)}, {'interface': IRGBColor, 'name': 'soil_color', 'value': (170, 85, 0)}, {'interface': 'IDict', 'name': 'colors', 'desc': 'dict (vid, rgb color) '}, {'interface': 'IBool', 'name': 'ambient_only', 'desc': 'If True, set to 0 all optical properties except the ambient one', 'value': False}],
                outputs=[{'interface': IInterface, 'name': 'scene', 'desc': 'PlantGL scene'}],
                widgetmodule=None,
                widgetclass=None,
               )

io_save_image = Factory(name='save_image',
                authors='C. Pradal, J. Coste ',
                description='Save a PlantGL scene in an image',
                category='io',
                nodemodule='image',
                nodeclass='save_image',
                inputs=[{'name': 'scene', 'desc': 'PlantGL Scene'}, 
                        {'interface': 'IStr', 'name': 'image_name', 'value': '%s/img%d.%s'}, 
                        {'interface': 'IDirStr', 'name': 'directory', 'value': '.'}, 
                        {'interface': 'IInt', 'name': 'index', 'value': 0}, 
                        {'interface': IEnumStr(enum=["png","jpg","tif"]), 'name': 'ext', 'value':'png'}],
                outputs=[{'interface': IInterface, 'name': 'scene', 'desc': 'PlantGL scene'}],
               )
__all__.append('io_save_image')





