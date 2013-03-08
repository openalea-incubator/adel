
# This file has been generated at Wed Mar  6 12:01:34 2013

from openalea.core import *


__name__ = 'alinea.adel.io'

__editable__ = True
__description__ = 'Manage input and output for adel such as meseasurement.'
__license__ = 'CECILL'
__url__ = ''
__alias__ = ['adel.io']
__version__ = '0.0.1'
__authors__ = 'C. Pradal, C. Fournier, C. Chambon'
__institutes__ = 'INRA, CIRAD, INRIA'
__icon__ = ''


__all__ = ['io_canL2canS', 'io_pandasDataframe2csv', 'io_select_adel_geometric_data', 'io_select_adel_botanic_data', 'io_dataframe', 'io_readRData', 'data2PandasDataframe', 'GetAdelString_GetAdelString', 'io_to_canestra', 'PairAsDict_PairAsDict', 'io_lpy2mtg', 'io_duplicate', 'io_mtg_factory', 'io_csvAsDict', 'io_saveRData', 'io_mtg2lpy', 'io_load_leaf_data', 'io_csv2pandasDataframe', 'io_dataframeAsdict', 'io_RlistAsDict', 'io_thermal_time', 'io_apply_property', 'image_save_image', 'pandasDataframe2data', 'io_to_plantgl']



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




io_pandasDataframe2csv = Factory(name='pandasDataframe2csv',
                authors='C. Chambon',
                description='Write a DataFrame to a comma-separated values (csv) file.',
                category='data i/o',
                nodemodule='io',
                nodeclass='pandasDataframe2csv',
                inputs=[{'interface': None, 'name': 'dataframe', 'value': None, 'desc': 'The DataFrame to write.'}, {'interface': IFileStr, 'name': 'csv_filepath', 'value': None, 'desc': 'The file path where the Dataframe is written.'}, {'interface': IStr, 'name': 'na_rep', 'value': None, 'desc': 'Missing data replacement.'}, {'interface': IBool, 'name': 'index', 'value': True, 'desc': 'Write row names (index)'}, {'interface': ISequence, 'name': 'index_label', 'value': None, 'desc': 'Column label for index column(s) if desired. If None is given, and header and index are True, then the index names are used. A sequence should be given if the DataFrame uses MultiIndex.'}],
                outputs=[{'interface': IFileStr, 'name': 'csv_filepath', 'desc': 'The file path where the Dataframe is written.'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_select_adel_geometric_data = Factory(name='select_adel_geometric_data',
                authors='C.Chambon',
                description='',
                category='data i/o',
                nodemodule='io',
                nodeclass='select_data_file',
                inputs=[{'interface': IDirStr, 'name': 'directory_path', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'laminaCur_pattern', 'value': 'laminaCur*.RData', 'desc': ''}, {'interface': IStr, 'name': 'lamina2D_pattern', 'value': 'lamina2D*.RData', 'desc': ''}, {'interface': IStr, 'name': 'laminaCur_filename', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'lamina2D_filename', 'value': None, 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'laminaCur_filepath'}, {'interface': IStr, 'name': 'lamina2D_filepath'}],
                widgetmodule='io',
                widgetclass='SelectDataFile',
               )




io_select_adel_botanic_data = Factory(name='select_adel_botanic_data',
                authors='C.Chambon',
                description='',
                category='data i/o',
                nodemodule='io',
                nodeclass='select_data_file',
                inputs=[{'interface': IDirStr, 'name': 'directory_path', 'value': None, 'desc': ''}, 
                        {'interface': IStr, 'name': 'axeT_pattern', 'value': '*axeT*.csv', 'desc': ''}, 
                        {'interface': IStr, 'name': 'dimT_pattern', 'value': '*dimT*.csv', 'desc': ''},
                        {'interface': IStr, 'name': 'phenT_pattern', 'value': '*phenT*.csv', 'desc': ''},
                        {'interface': IStr, 'name': 'earT_pattern', 'value': '*earT*.csv', 'desc': ''},
                        {'interface': IStr, 'name': 'ssi2sen_pattern', 'value': '*ssi2sen*.csv', 'desc': ''},
                        {'interface': IStr, 'name': 'axeT_filename', 'value': None, 'desc': ''}, 
                        {'interface': IStr, 'name': 'dimT_filename', 'value': None, 'desc': ''},
                        {'interface': IStr, 'name': 'phenT_filename', 'value': None, 'desc': ''},
                        {'interface': IStr, 'name': 'earT_filename', 'value': None, 'desc': ''},
                        {'interface': IStr, 'name': 'ssi2sen_filename', 'value': None, 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'axeT_filepath'}, 
                         {'interface': IStr, 'name': 'dimT_filepath'},
                         {'interface': IStr, 'name': 'phenT_filepath'},
                         {'interface': IStr, 'name': 'earT_filepath'},
                         {'interface': IStr, 'name': 'ssi2sen_filepath'}],
                widgetmodule='io',
                widgetclass='SelectDataFile',
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




data2PandasDataframe = CompositeNodeFactory(name='data2PandasDataframe',
                             description='Composite node made by grouping openalea.data file.get_data and alinea.adel.io.csv2pandasDataframe',
                             category='data i/o',
                             doc='',
                             inputs=[  {  'desc': '', 'interface': IStr, 'name': 'package(get_data)', 'value': None},
   {  'desc': '', 'interface': IStr, 'name': 'glob(get_data)', 'value': '*'},
   {  'desc': '',
      'interface': IStr,
      'name': 'filename(get_data)',
      'value': None}],
                             outputs=[  {  'desc': 'A pandas.DataFrame instance which represents the csv file.',
      'interface': None,
      'name': 'dataframe(csv2pandasDataframe)'}],
                             elt_factory={  19: ('openalea.data file', 'get_data'),
   20: ('alinea.adel.io', 'csv2pandasDataframe')},
                             elt_connections={  34301880: ('__in__', 2, 19, 2),
   34301904: ('__in__', 1, 19, 1),
   34301928: ('__in__', 0, 19, 0),
   34301952: (20, 0, '__out__', 0),
   34302048: (19, 0, 20, 0)},
                             elt_data={  19: {  'block': False,
          'caption': 'get_data',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x3dc71d0> : "get_data"',
          'hide': True,
          'id': 19,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 509.64542300600647,
          'posy': -31.791839043245815,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   20: {  'block': False,
          'caption': 'csv2pandasDataframe',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x3db4910> : "csv2pandasDataframe"',
          'hide': True,
          'id': 20,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 491.8627065927451,
          'posy': 21.13472567390702,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'delay': 0,
                'hide': True,
                'id': 0,
                'lazy': True,
                'port_hide_changed': set(),
                'posx': 521.8488883496688,
                'posy': -88.5714957876693,
                'priority': 0,
                'use_user_color': False,
                'user_application': None,
                'user_color': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'delay': 0,
                 'hide': True,
                 'id': 1,
                 'lazy': True,
                 'port_hide_changed': set(),
                 'posx': 550.9668572660541,
                 'posy': 78.82221022550736,
                 'priority': 0,
                 'use_user_color': False,
                 'user_application': None,
                 'user_color': None}},
                             elt_value={  19: [],
   20: [(1, 'None'), (2, 'None'), (3, 'False')],
   '__in__': [],
   '__out__': []},
                             elt_ad_hoc={  19: {  'position': [509.64542300600647, -31.791839043245815],
          'useUserColor': False,
          'userColor': None},
   20: {  'position': [491.8627065927451, 21.13472567390702],
          'useUserColor': False,
          'userColor': None},
   '__in__': {  'position': [521.8488883496688, -88.5714957876693],
                'useUserColor': False,
                'userColor': None},
   '__out__': {  'position': [550.9668572660541, 78.82221022550736],
                 'useUserColor': False,
                 'userColor': None}},
                             lazy=True,
                             eval_algo='LambdaEvaluation',
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




io_csv2pandasDataframe = Factory(name='csv2pandasDataframe',
                authors='C. Chambon',
                description='Read CSV (comma-separated) file into DataFrame.',
                category='data i/o',
                nodemodule='io',
                nodeclass='csv2pandasDataframe',
                inputs=[{'interface': IFileStr, 'name': 'csv_filepath', 'value': None, 'desc': 'The filepath of the csv file to import.'}, {'interface': ISequence, 'name': 'index_col', 'value': None, 'desc': 'Column to use as the row labels of the DataFrame. If a sequence is given, a MultiIndex is used.'}, {'interface': ISequence, 'name': 'na_values', 'value': None, 'desc': 'List of additional strings to recognize as NA/NaN.'}, {'interface': IBool, 'name': 'parse_dates', 'value': False, 'desc': 'Attempt to parse dates in the index column(s).'}],
                outputs=[{'interface': None, 'name': 'dataframe', 'desc': 'A pandas.DataFrame instance which represents the csv file.'}],
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




image_save_image = Factory(name='save_image',
                authors='C. Pradal, J. Coste ',
                description='Save a PlantGL scene in an image',
                category='io',
                nodemodule='image',
                nodeclass='save_image',
                inputs=[{'name': 'scene', 'desc': 'PlantGL Scene'}, {'interface': 'IStr', 'name': 'image_name', 'value': '%s/img%04d.%s'}, {'interface': 'IDirStr', 'name': 'directory', 'value': '.'}, {'interface': 'IInt', 'name': 'index', 'value': 0}, {'interface': IEnumStr(enum=['png', 'jpg', 'tif']), 'name': 'ext', 'value': 'png'}],
                outputs=[{'interface': IInterface, 'name': 'scene', 'desc': 'PlantGL scene'}],
                widgetmodule=None,
                widgetclass=None,
               )




pandasDataframe2data = CompositeNodeFactory(name='pandasDataframe2data',
                             description='Composite node made by grouping openalea.data file.get_data and alinea.adel.io.pandasDataframe2csv',
                             category='data i/o',
                             doc='',
                             inputs=[  {  'desc': 'The DataFrame to write.',
      'interface': None,
      'name': 'dataframe(pandasDataframe2csv)',
      'value': None},
   {  'desc': '', 'interface': IStr, 'name': 'package(get_data)', 'value': None},
   {  'desc': '', 'interface': IStr, 'name': 'glob(get_data)', 'value': '*'},
   {  'desc': '',
      'interface': IStr,
      'name': 'filename(get_data)',
      'value': None}],
                             outputs=[  {  'desc': 'The file path where the Dataframe is written.',
      'interface': IFileStr,
      'name': 'csv_filepath(pandasDataframe2csv)'}],
                             elt_factory={  8: ('openalea.data file', 'get_data'),
   14: ('alinea.adel.io', 'pandasDataframe2csv')},
                             elt_connections={  20768744: ('__in__', 2, 8, 1),
   20768768: (14, 0, '__out__', 0),
   20768792: ('__in__', 3, 8, 2),
   20768816: (8, 0, 14, 1),
   20768840: ('__in__', 0, 14, 0),
   20768864: ('__in__', 1, 8, 0)},
                             elt_data={  8: {  'block': False,
         'caption': 'get_data',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x30c8d50> : "get_data"',
         'hide': True,
         'id': 8,
         'lazy': True,
         'port_hide_changed': set(),
         'posx': 371.164776354299,
         'posy': -25.54821877276729,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   14: {  'block': False,
          'caption': 'pandasDataframe2csv',
          'delay': 0,
          'factory': '<openalea.core.node.NodeFactory object at 0x30c8410> : "pandasDataframe2csv"',
          'hide': True,
          'id': 14,
          'lazy': True,
          'port_hide_changed': set(),
          'posx': 305.4928529228377,
          'posy': 48.03346744069775,
          'priority': 0,
          'use_user_color': False,
          'user_application': None,
          'user_color': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'delay': 0,
                'hide': True,
                'id': 0,
                'lazy': True,
                'port_hide_changed': set(),
                'posx': 370.64033370419804,
                'posy': -204.7760830037442,
                'priority': 0,
                'use_user_color': False,
                'user_application': None,
                'user_color': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'delay': 0,
                 'hide': True,
                 'id': 1,
                 'lazy': True,
                 'port_hide_changed': set(),
                 'posx': 396.3969408123752,
                 'posy': 168.44012208426494,
                 'priority': 0,
                 'use_user_color': False,
                 'user_application': None,
                 'user_color': None}},
                             elt_value={  8: [],
   14: [(2, "'NA'"), (3, 'False'), (4, 'None')],
   '__in__': [],
   '__out__': []},
                             elt_ad_hoc={  8: {'position': [371.164776354299, -25.54821877276729], 'userColor': None, 'useUserColor': False},
   14: {'position': [305.4928529228377, 48.03346744069775], 'userColor': None, 'useUserColor': False},
   '__in__': {'position': [370.64033370419804, -204.7760830037442], 'userColor': None, 'useUserColor': False},
   '__out__': {'position': [396.3969408123752, 168.44012208426494], 'userColor': None, 'useUserColor': False}},
                             lazy=True,
                             eval_algo='LambdaEvaluation',
                             )




io_to_plantgl = Factory(name='to_plantgl',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Adapt a Canestra scene to a PlantGL scene',
                category='io',
                nodemodule='io',
                nodeclass='to_plantgl',
                inputs=[{'name': 'scene', 'desc': 'Canestra Scene'}, {'interface': IRGBColor, 'name': 'leaf_color', 'value': (0, 180, 0)}, {'interface': IRGBColor, 'name': 'stem_color', 'value': (0, 130, 0)}, {'interface': IRGBColor, 'name': 'soil_color', 'value': (170, 85, 0)}, {'interface': 'IDict', 'name': 'colors', 'desc': 'dict (vid, rgb color) '}, {'interface': 'IBool', 'name': 'ambient_only', 'value': False, 'desc': 'If True, set to 0 all optical properties except the ambient one'}],
                outputs=[{'interface': IInterface, 'name': 'scene', 'desc': 'PlantGL scene'}],
                widgetmodule=None,
                widgetclass=None,
               )




