# This file stores a set of data which can be used by the function 
# alinea.adel.plantgen.read_plantgen_inputs to define the inputs of the function 
# alinea.adel.plantgen.gen_adel_input_data.
# dynT_user and dimT_user are absolute file paths.
# ATTENTION: this file is just an example. Please adapt manually the 
# value of dynT_user and dimT_user to your need.
dynT_user = 'c://openaleapkg//adel//adel//data//Maxwell_2011_d220_N1_dynT_user_SHORT.csv'
dimT_user = 'c://openaleapkg//adel//adel//data//Maxwell_2011_d220_N1_dimT_user_SHORT.csv'
plants_number = 217
plants_density = 217
decide_child_axis_probabilities = {'T0': 0.0, 'T1': 0.883333333333333,
                                   'T2': 0.966666666666667, 'T3': 0.816666666666667,
                                   'T4': 0.0833333333333333}
MS_leaves_number_probabilities = {'10': 0.146, 
                                  '11': 0.818,
                                  '12': 0.036,
                                  '13': 0.0,
                                  '14': 0.0}
ears_density = 344
GL_number = {1117.03894252306: 5.56702696865511, 1212.13219252307:5.31217286737119,
             1368.69019252307:4.75030723720174, 1686.7860258564:2.52390175114029 , 1818.0:0.0 }
delais_TT_stop_del_axis = 600
TT_col_break = 0.0
inner_params = {'DELAIS_PHYLL_COL_TIP_1ST': 1.0,
                'DELAIS_PHYLL_COL_TIP_NTH': 1.6}