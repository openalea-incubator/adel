# This file stores a set of data which can be used by the function 
# alinea.adel.plantgen.read_plantgen_inputs to define the inputs of the function 
# alinea.adel.plantgen.gen_adel_input_data_from_min.
# ATTENTION: this file is just an example and IS NOT portable. Please adapt it to 
# your need manually or use the graphic node read_plantgen_inputs instead.
dynT_user = {'a_cohort': 0.0102,
             'TT_col_0': -0.771289027,
             'n0': 4.871559739,
             'n1': 3.24283148,
             'n2': 5.8,
             'TT_col_N_phytomer_potential': {'MS': 1078.0,
                                             'T1': 1148.0,
                                             'T2': 1158.0,
                                             'T3': 1168.0,
                                             'T4': 1178.0}}
dimT_user = '/home/cchambon/workspace/openaleapkg_tr/adel/test/data/test_plantgen/MIN_MIN/dimT_user.csv'
plants_number = 100
plants_density = 250
decide_child_axis_probabilities = {'T0': 0.0, 'T1': 0.900,
                                   'T2': 0.983, 'T3': 0.817,
                                   'T4': 0.117}
MS_leaves_number_probabilities = {'10': 0.145,
                                  '11': 0.818,
                                  '12': 0.036,
                                  '13': 0.0,
                                  '14': 0.0}
ears_density = 500
GL_number = {1117.0: 5.6, 1212.1:5.4,
             1368.7:4.9, 1686.8:2.4,
             1880.0:0.0}
delais_TT_stop_del_axis = 600
TT_col_break = 0.0