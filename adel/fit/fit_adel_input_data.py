from adel.fit.fit_adel_input_data_first import fit_adel_input_data_first
from adel.fit.fit_adel_input_data_second import fit_adel_input_data_second

class DataCompleteness:
    MIN=1
    SHORT=2
    FULL=3


def fit_adel_input_data(user_parameters,
                        user_dims,
                        plant_number=100, 
                        cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.983, '6': 0.817, '7': 0.117, '8': 0.0, '9': 0.0, '10': 0.0}, 
                        main_stem_leaves_number_probability_distribution={'10': 0.145, '11': 0.818, '12': 0.036, '13': 0.0, '14': 0.0},
                        bolting_date=500,
                        flowering_date=1440,
                        final_axes_number=250,
                        GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                        delais_TT_stop_del_axis=600,
                        dTT_MS_cohort={'4': 70, '5': 80, '6': 90, '7': 100},
                        TT_col_break=0.0,
                        user_parameters_completeness=DataCompleteness.MIN,
                        user_dims_completeness=DataCompleteness.MIN,
                        ):
    '''
    Fit ADEL input data.
    
    :Parameters:
        - `user_parameters` : the parameters set by the user. Can be either a dict or a table.  
        - `user_dims` : a table describing the organ dimensions set by the user.
        - `plant_number` : the number of plants.
        - `cohort_probabilities` : the probability of emergence of a child axis when the parent axis is present. This probability is 
           related to the cohort of the child axis.  
        - `main_stem_leaves_number_probability_distribution` : the probability distribution of the main stem leaves number.
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
        - `final_axes_number` : the final number of axes which have an ear, per square meter.
        - `GL_number` : the GL decimal number measured at several thermal time (including the senescence end).
        - `delais_TT_stop_del_axis` : Thermal time during which a tiller remains present on the plant after the tiller has stopped growing.
        - `dTT_MS_cohort` : ???
        - `TT_col_break` : ???
        - `user_parameters_completeness` : the level of completeness of the parameters set by the user. See adel.fit.fit_adel_input_data.DataCompleteness.
        - `user_dims_completeness` : the level of completeness of the organ dimensions set by the user. See adel.fit.fit_adel_input_data.DataCompleteness.
        
    :Types:
        - `user_parameters` : dict | pandas.DataFrame      
        - `user_dims` : pandas.DataFrame
        - `plant_number` : int
        - `cohort_probabilities` : dict of str:float
        - `main_stem_leaves_number_probability_distribution` : dict of str:float
        - `bolting_date` : int
        - `flowering_date` : int
        - `final_axes_number` : int
        - `GL_number` : dict of float:float
        - `delais_TT_stop_del_axis` : int
        - `dTT_MS_cohort` : dict of str:float
        - `TT_col_break` : float 
        - `user_parameters_completeness` : int
        - `user_dims_completeness` : int
        
    :return: the axis table, the absolute phen table, the relative phen table, \
           the absolute dim table, the parameters table, the first leaf phen table, \
           the HS_GL_SSI dynamic table, the relative dim table, and \
           the tillering dynamic table
    
    :rtype: tuple of pandas.DataFrame
    '''
        
    (first_axis_table_dataframe, 
    dim_table_dataframe, 
    parameters_table_dataframe, 
    tillering_dynamic_dataframe) = fit_adel_input_data_first(plant_number, 
                                                             cohort_probabilities, 
                                                             main_stem_leaves_number_probability_distribution, 
                                                             bolting_date, 
                                                             flowering_date, 
                                                             final_axes_number) 
    parameters_table_dataframe.ix[0]['TT_col_break'] = TT_col_break   
    if user_parameters_completeness == DataCompleteness.MIN:
        parameters_table_dataframe.ix[0]['a_cohort'] = user_parameters['a_cohort']
        parameters_table_dataframe.ix[0]['TT_col_0'] = user_parameters['TT_col_0']
        parameters_table_dataframe.ix[0]['TT_col_nff'] = user_parameters['TT_col_nff']
        for N_cohort, parameters_table_dataframe_grouped_by_N_cohort in parameters_table_dataframe.groupby('N_cohort'):
            N_cohort_int = int(N_cohort)
            if N_cohort_int == 1:
                current_dTT_MS_cohort = 0.0
            else:
                current_dTT_MS_cohort = dTT_MS_cohort[str(N_cohort_int)]
            index_to_set = parameters_table_dataframe_grouped_by_N_cohort.index[0]
            parameters_table_dataframe['dTT_MS_cohort'][index_to_set] = current_dTT_MS_cohort        
        parameters_table_dataframe.ix[0]['n0'] = user_parameters['n0']
        parameters_table_dataframe.ix[0]['n1'] = user_parameters['n1']
        parameters_table_dataframe.ix[0]['n2'] = user_parameters['n2']
    elif user_parameters_completeness == DataCompleteness.SHORT:
        user_grouped = user_parameters.groupby('N_cohort')
        for N_cohort, generated_group in parameters_table_dataframe.groupby('N_cohort'):
            if N_cohort not in user_grouped.groups:
                continue
            index_to_get = user_grouped.get_group(N_cohort).index[0]
            index_to_set = generated_group.index[0]
            columns_to_set = user_parameters.columns
            parameters_table_dataframe.ix[index_to_set][columns_to_set] = user_parameters.ix[index_to_get]        
    elif user_parameters_completeness == DataCompleteness.FULL:
        user_grouped = user_parameters.groupby(['N_cohort', 'Nff'])
        for (N_cohort, Nff), generated_group in parameters_table_dataframe.groupby(['N_cohort', 'Nff']):
            if (N_cohort, Nff) not in user_grouped.groups:
                continue
            index_to_get = user_grouped.get_group((N_cohort, Nff)).index[0]
            index_to_set = generated_group.index[0]
            columns_to_set = user_parameters.columns
            parameters_table_dataframe.ix[index_to_set][columns_to_set] = user_parameters.ix[index_to_get]
    else:
        raise Exception('''%s is an unknown user data completeness value. \
The values can be one of %s''', (str(user_parameters_completeness), 
                                 str([DataCompleteness.MIN, 
                                      DataCompleteness.SHORT, 
                                      DataCompleteness.FULL])))
            
    organ_dim_list = ['L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
    if user_dims_completeness == DataCompleteness.MIN:
        index_to_set = dim_table_dataframe[dim_table_dataframe['id_dim']==dim_table_dataframe['id_dim'][0]].index
        for organ_dim in organ_dim_list:
            dim_table_dataframe[organ_dim][index_to_set] = user_dims[organ_dim][index_to_set]
    elif user_dims_completeness == DataCompleteness.SHORT:
        user_grouped = user_dims.groupby('id_axis')
        for generated_id_dim, generated_group in dim_table_dataframe.groupby('id_dim'):
            generated_id_axis = float(str(int(generated_id_dim))[:-2])
            if generated_id_axis not in user_grouped.groups:
                continue
            user_group = user_grouped.get_group(generated_id_axis)
            user_index_phytomer = user_group['index_phytomer'][user_group.index[-1]]
            user_id_dim = float(''.join([str(int(generated_id_axis)), str(int(user_index_phytomer))]))
            if user_id_dim != generated_id_dim:
                continue
            index_to_get = user_group.index
            index_to_set = generated_group.index
            for organ_dim in organ_dim_list:
                dim_table_dataframe[organ_dim][index_to_set] = user_dims[organ_dim][index_to_get]     
    elif user_dims_completeness == DataCompleteness.FULL:
        user_grouped = user_dims.groupby('id_dim')
        for generated_id_dim, generated_group in dim_table_dataframe.groupby('id_dim'):
            if generated_id_dim not in user_grouped.groups:
                continue
            user_group = user_grouped.get_group(generated_id_dim)
            dim_table_dataframe.ix[generated_group.index] = user_dims.ix[user_group.index]
    else:
        raise Exception('''%s is an unknown user data completeness value. \
The values can be one of %s''', (str(user_dims_completeness), 
                                 str([DataCompleteness.MIN, 
                                      DataCompleteness.SHORT, 
                                      DataCompleteness.FULL])))
        
    (second_axis_table_dataframe, 
     absolute_second_phen_table_dataframe, 
     relative_second_phen_table_dataframe,
     absolute_dim_table_dataframe, 
     second_parameters_table_dataframe, 
     first_leaf_phen_table_dataframe,
     HS_GL_SSI_dynamic_dataframe, 
     relative_dim_table_dataframe) = fit_adel_input_data_second(first_axis_table_dataframe, 
                                                                dim_table_dataframe, 
                                                                parameters_table_dataframe, 
                                                                GL_number, 
                                                                bolting_date, 
                                                                flowering_date, 
                                                                delais_TT_stop_del_axis, 
                                                                final_axes_number)
    
    return second_axis_table_dataframe, absolute_second_phen_table_dataframe, relative_second_phen_table_dataframe, \
           absolute_dim_table_dataframe, second_parameters_table_dataframe, first_leaf_phen_table_dataframe, \
           HS_GL_SSI_dynamic_dataframe, relative_dim_table_dataframe, tillering_dynamic_dataframe
    
    

