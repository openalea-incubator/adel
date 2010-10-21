
# This file has been generated at Wed Aug 04 16:49:56 2010

from openalea.core import *


__name__ = 'alinea.nema.modules'

__editable__ = True
__description__ = ''
__license__ = ''
__url__ = ''
__alias__ = []
__version__ = ''
__authors__ = ''
__institutes__ = ''
__icon__ = ''


__all__ = ['sheath_sheath', 'internode_internode', 'root_root', 'peduncle_peduncle', 'grain_grain', 'chaff_chaff', 'lamina_lamina']



sheath_sheath = Factory(name='sheath',
                description='create a sheath object for Nema',
                category='model',
                nodemodule='sheath',
                nodeclass='sheath',
                inputs=[{'interface': IDict, 'name': 'Parameters', 'value': {'pm1': 0.80000000000000004, 'pm2': 6.0, 'peff': 5.0000000000000004e-06, 'pdeath': 0.40000000000000002, 'beta': 2.0, 'SynthRate': 0.00014999999999999999, 'alpha': 2.0, 'DegRate': 0.0080000000000000002, 'p': 0.10000000000000001, 'InsertionAngle': 1.5700000000000001, 'k2': 864000.0, 'k1': 0.0018, 'DegDMRate': 0.0080000000000000002, 'TTexp': 940.0, 'sink': 1.0}, 'desc': ''}, {'interface': IDict, 'name': 'Initial state', 'value': {'Nphs': [2.5999999999999998e-05, 0.00018000000000000001, 0.00066, 0.0018], 'Area': [0.00020000000000000001, 0.00040000000000000002, 0.00050000000000000001, 0.00059999999999999995], 'TTinit': [-510.0, -420.0, -330.0, -240.0], 'DMrem': [0.0019, 0.0092999999999999992, 0.014800000000000001, 0.022200000000000001], 'Length': [0.11, 0.125, 0.14000000000000001, 0.14499999999999999], 'Nstruct': [0.00012, 0.00011, 0.00021000000000000001, 0.00068000000000000005], 'DMstruct': [0.0086, 0.042799999999999998, 0.068500000000000005, 0.1027]}, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'sheath object', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




internode_internode = Factory(name='internode',
                description='create a internode object for Nema',
                category='model',
                nodemodule='internode',
                nodeclass='internode',
                inputs=[{'interface': IDict, 'name': 'Parameters', 'value': {'pm1': 0.80000000000000004, 'pm2': 6.0, 'peff': 5.0000000000000004e-06, 'pdeath': 0.40000000000000002, 'beta': 2.0, 'SynthRate': 3.0000000000000001e-05, 'alpha': 2.0, 'DegRate': 0.0080000000000000002, 'p': 0.10000000000000001, 'InsertionAngle': 1.5700000000000001, 'k2': 864000.0, 'k1': 0.0018, 'DegDMRate': 0.0080000000000000002, 'TTexp': 920.0, 'sink': 1.0}, 'desc': ''}, {'interface': IDict, 'name': 'Initial state', 'value': {'Nphs': [0.0, 0.0, 0.00027999999999999998, 0.00123], 'Area': [0.0001, 0.00025000000000000001, 0.00040000000000000002, 0.0015], 'TTinit': [-490.0, -400.0, -310.0, -220.0], 'DMrem': [0.0092999999999999992, 0.033000000000000002, 0.039, 0.041000000000000002], 'Length': [0.050000000000000003, 0.085999999999999993, 0.128, 0.186], 'Nstruct': [5.0000000000000002e-05, 0.00013999999999999999, 0.00033, 0.00084999999999999995], 'DMstruct': [0.042799999999999998, 0.15409999999999999, 0.1797, 0.188]}, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'internode object', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




root_root = Factory(name='root',
                description='create a root object for Nema',
                category='model',
                nodemodule='root',
                nodeclass='root',
                inputs=[{'interface': IDict, 'name': 'Parameters', 'value': {'TTinit': -1000.0, 'kroot2': 5.0000000000000004e-06, 'kroot1': 2.5, 'beta': 2.0, 'sink': 1.0, 'QonDmin': 0.0012999999999999999, 'alpha': 2.0, 'degDMcoef': 0.0080000000000000002, 'NsoilMin': 0.0, 'degcoef': 0.0080000000000000002, 'p': 0.10000000000000001, 'beta2': 2300.0, 'beta1': 130.0, 'TTexp': 1700.0, 'Umax': 0.0001}, 'desc': ''}, {'interface': IDict, 'name': 'Initial state', 'value': {'DMrem': 0.056000000000000001, 'Nstruct': 0.01, 'Nrem': 0.0011000000000000001, 'DMstruct': 0.504}, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'root object', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




peduncle_peduncle = Factory(name='peduncle',
                description='create a peduncle object for Nema',
                category='model',
                nodemodule='peduncle',
                nodeclass='peduncle',
                inputs=[{'interface': IDict, 'name': 'Parameters', 'value': {'pm1': 0.80000000000000004, 'pm2': 6.0, 'peff': 5.0000000000000004e-06, 'pdeath': 0.40000000000000002, 'beta': 2.0, 'SynthRate': 3.0000000000000001e-05, 'alpha': 2.0, 'DegRate': 0.0080000000000000002, 'p': 0.10000000000000001, 'InsertionAngle': 1.5700000000000001, 'k2': 864000.0, 'k1': 0.0018, 'DegDMRate': 0.0080000000000000002, 'TTexp': 830.0, 'sink': 2.0}, 'desc': ''}, {'interface': IDict, 'name': 'Initial state', 'value': {'Nphs': [0.0025300000000000001], 'Area': [0.0023999999999999998], 'TTinit': [-130.0], 'DMrem': [0.055500000000000001], 'Length': [0.219], 'Nstruct': [0.0012999999999999999], 'DMstruct': [0.25679999999999997]}, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'peduncle object', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




grain_grain = Factory(name='grain',
                description='create a grain object for Nema',
                category='model',
                nodemodule='grain',
                nodeclass='grain',
                inputs=[{'interface': IDict, 'name': 'Parameters', 'value': {'p': 6.0, 'beta': 2.0, 'sink': 7.0, 'alpha': 2.0, 'TTinit': 0.0, 'TTexp': 700.0, 'tgrain': 250.0, 'ggrain': 0.0054999999999999997}, 'desc': ''}, {'interface': IDict, 'name': 'Initial state', 'value': {'Ngrain': 0.0023999999999999998, 'DMgrain': 0.12}, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'grain object', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




chaff_chaff = Factory(name='chaff',
                description='create a chaff object for Nema',
                category='model',
                nodemodule='chaff',
                nodeclass='chaff',
                inputs=[{'interface': IDict, 'name': 'Parameters', 'value': {'pm1': 0.80000000000000004, 'pm2': 6.0, 'peff': 5.0000000000000004e-06, 'pdeath': 0.40000000000000002, 'beta': 2.0, 'SynthRate': 3.0000000000000001e-05, 'alpha': 2.0, 'DegRate': 0.0080000000000000002, 'p': 0.10000000000000001, 'InsertionAngle': 1.5700000000000001, 'k2': 864000, 'k1': 0.0018, 'DegDMRate': 0.0080000000000000002, 'TTexp': 1100.0, 'sink': 2.0}, 'desc': ''}, {'interface': IDict, 'name': 'Initial state', 'value': {'Nphs': [0.0036800000000000001], 'Area': [0.00075000000000000002], 'TTinit': [-400], 'DMrem': [0.050000000000000003], 'Length': [0.089999999999999997], 'Nstruct': [0.00107], 'DMstruct': [0.20999999999999999]}, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'chaff object', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




lamina_lamina = Factory(name='lamina',
                description='create a lamina object for Nema',
                category='model',
                nodemodule='lamina',
                nodeclass='lamina',
                inputs=[{'interface': IDict, 'name': 'Parameters', 'value': {'pm1': 0.80000000000000004, 'pm2': 6.0, 'peff': 5.0000000000000004e-06, 'pdeath': 0.40000000000000002, 'beta': 2.0, 'SynthRate': 0.00014999999999999999, 'alpha': 2.0, 'DegRate': 0.0080000000000000002, 'p': 0.10000000000000001, 'InsertionAngle': 0.52000000000000002, 'k2': 864000.0, 'k1': 0.0018, 'DegDMRate': 0.0080000000000000002, 'TTexp': 1090.0, 'sink': 1.0}, 'desc': ''}, {'interface': None, 'name': 'Initial state', 'value': {'Nphs': [9.0000000000000006e-05, 0.0011999999999999999, 0.0028999999999999998, 0.0053], 'Area': [0.0016000000000000001, 0.0022799999999999999, 0.0033999999999999998, 0.00346], 'TTinit': [-660.0, -570.0, -480.0, -390.0], 'DMrem': [0.0, 0.02, 0.040000000000000001, 0.040000000000000001], 'Length': [0.182, 0.21099999999999999, 0.22650000000000001, 0.17349999999999999], 'Nstruct': [0.00038000000000000002, 0.00052999999999999998, 0.00083000000000000001, 0.0010200000000000001], 'DMstruct': [0.040000000000000001, 0.050000000000000003, 0.089999999999999997, 0.14000000000000001]}, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'lamina object', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




