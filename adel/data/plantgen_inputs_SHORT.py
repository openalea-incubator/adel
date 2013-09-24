# This file stores a set of data which can be used by the function 
# alinea.adel.plantgen.read_plantgen_inputs to define the inputs of the function 
# alinea.adel.plantgen.gen_adel_input_data.
# dynT_user and dimT_user are absolute file paths.
# ATTENTION: this file is just an example. Please adapt manually the 
# value of dynT_user and dimT_user to your need.
dynT_user = '/home/cchambon/workspace/openaleapkg_tr/adel/adel/data/Mariem_dynT_user_SHORT.csv'
dimT_user = '/home/cchambon/workspace/openaleapkg_tr/adel/adel/data/Mariem_dimT_user_SHORT.csv'
plants_number = 100
plants_density = 250
decide_child_axis_probabilities = {'T0': 0.0, 'T1': 0.900,
                                   'T2': 0.983, 'T3': 0.817,
                                   'T4': 0.117}
MS_leaves_number_probabilities = {'10': 0.145,
                                  '11': 0.818,
                                  '12': 0.037,
                                  '13': 0.0,
                                  '14': 0.0}
ears_density = 500
GL_number = {1117.0: 5.6, 1212.1:5.4,
             1368.7:4.9, 1686.8:2.4,
             1880.0:0.0}
delais_TT_stop_del_axis = 600
TT_col_break = 0.0
inner_params = {'DELAIS_PHYLL_COL_TIP_1ST': 1.0,
                'DELAIS_PHYLL_COL_TIP_NTH': 1.6}