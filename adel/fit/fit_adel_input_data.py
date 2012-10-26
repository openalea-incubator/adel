# -*- python -*-
#
#       Adel.Fit
#
#       Copyright 2006-2012 INRIA - CIRAD - INRA
#
#       File author(s): Camille Chambon <camille.chambon@grignon.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################
import numpy as np
import pandas

from adel.fit.fit_adel_input_data_first import fit_adel_input_data_first
from adel.fit.fit_adel_input_data_second import fit_adel_input_data_second

class DataCompleteness:
    MIN=1
    SHORT=2
    FULL=3


def fit_adel_input_data(user_leaf_dynamic_parameters,
                        user_dims,
                        plant_number=100, 
                        cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.983, '6': 0.817, '7': 0.117}, 
                        main_stem_leaves_number_probability_distribution={'10': 0.145, '11': 0.818, '12': 0.036, '13': 0.0, '14': 0.0},
                        bolting_date=500,
                        flowering_date=1440,
                        final_axes_number=250,
                        GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                        delais_TT_stop_del_axis=600,
                        TT_col_break=0.0,
                        user_leaf_dynamic_parameters_completeness=DataCompleteness.MIN,
                        user_dims_completeness=DataCompleteness.MIN,
                        ):
    '''
    Fit ADEL input data.
    
    :Parameters:
        - `user_leaf_dynamic_parameters` : the leaf_dynamic_parameters set by the user. Can be either a dict or a table.  
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
        - `user_leaf_dynamic_parameters_completeness` : the level of completeness of the leaf_dynamic_parameters set by the user. See adel.fit.fit_adel_input_data.DataCompleteness.
        - `user_dims_completeness` : the level of completeness of the organ dimensions set by the user. See adel.fit.fit_adel_input_data.DataCompleteness.
        
    :Types:
        - `user_leaf_dynamic_parameters` : dict | pandas.DataFrame      
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
        - `user_leaf_dynamic_parameters_completeness` : int
        - `user_dims_completeness` : int
        
    :return: the axis table, the absolute phen table, the relative phen table, \
           the absolute dim table, the leaf_dynamic_parameters table, the first leaf phen table, \
           the HS_GL_SSI dynamic table, the relative dim table, and \
           the tillering dynamic table
    
    :rtype: tuple of pandas.DataFrame
    '''
    # 1. check input leaf_dynamic_parameters validity
    assert user_leaf_dynamic_parameters_completeness in DataCompleteness.__dict__.values()
    assert user_dims_completeness in DataCompleteness.__dict__.values()
    assert isinstance(TT_col_break, float)
    if user_leaf_dynamic_parameters_completeness == DataCompleteness.MIN:
        # check user_leaf_dynamic_parameters validity
        expected_user_leaf_dynamic_parameters_keys_value_types = {'a_cohort': float, 
                                                     'TT_col_0': float, 
                                                     'TT_col_nff': dict, 
                                                     'n0': float,
                                                     'n1': float,
                                                     'n2': float}
        user_leaf_dynamic_parameters_keys_value_types = dict(zip(user_leaf_dynamic_parameters.keys(), 
                                                    [type(value) for value in user_leaf_dynamic_parameters.values()]))
        assert expected_user_leaf_dynamic_parameters_keys_value_types == user_leaf_dynamic_parameters_keys_value_types
        # check user_dims validity
        assert isinstance(user_dims, pandas.DataFrame)
        expected_user_dims_columns = ['index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        assert (user_dims.columns == expected_user_dims_columns).all()
        assert user_dims.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        assert user_dims['index_phytomer'].unique().size == user_dims['index_phytomer'].size
    elif user_leaf_dynamic_parameters_completeness == DataCompleteness.SHORT:
        # check user_leaf_dynamic_parameters validity
        assert isinstance(user_leaf_dynamic_parameters, pandas.DataFrame)
        expected_user_leaf_dynamic_parameters_columns = ['N_cohort', 'a_cohort', 'TT_col_0', 'TT_col_nff', 'n0', 'n1', 'n2']
        assert (user_leaf_dynamic_parameters.columns == expected_user_leaf_dynamic_parameters_columns).all()
        assert user_leaf_dynamic_parameters.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        assert user_leaf_dynamic_parameters['N_cohort'].unique().size == user_leaf_dynamic_parameters['N_cohort'].size
        # check user_dims validity
        assert isinstance(user_dims, pandas.DataFrame)
        expected_user_dims_columns = ['id_axis', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        assert (user_dims.columns == expected_user_dims_columns).all()
        assert user_dims.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        grouped = user_dims.groupby(['id_axis', 'index_phytomer'])
        assert len(grouped.groups) == user_dims.index.size            
    elif user_leaf_dynamic_parameters_completeness == DataCompleteness.FULL:
        # check user_leaf_dynamic_parameters validity
        assert isinstance(user_leaf_dynamic_parameters, pandas.DataFrame)
        expected_user_leaf_dynamic_parameters_columns = ['N_cohort', 'Nff', 'a_cohort', 'TT_col_0', 'TT_col_nff', 'n0', 'n1', 'n2']
        assert (user_leaf_dynamic_parameters.columns == expected_user_leaf_dynamic_parameters_columns).all()
        assert user_leaf_dynamic_parameters.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        grouped = user_leaf_dynamic_parameters.groupby(['N_cohort', 'Nff'])
        assert len(grouped.groups) == user_leaf_dynamic_parameters.index.size    
        # check user_dims validity
        assert isinstance(user_dims, pandas.DataFrame)
        expected_user_dims_columns = ['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        assert (user_dims.columns == expected_user_dims_columns).all()
        assert user_dims.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        grouped = user_dims.groupby(['id_dim', 'index_phytomer'])
        assert len(grouped.groups) == user_dims.index.size    
    assert isinstance(plant_number, int)
    assert isinstance(cohort_probabilities, dict)
    assert isinstance(main_stem_leaves_number_probability_distribution, dict)
    assert isinstance(bolting_date, int)
    assert isinstance(flowering_date, int)
    assert isinstance(final_axes_number, int)
    assert isinstance(GL_number, dict)
    assert isinstance(delais_TT_stop_del_axis, int) 
    
    # 2. first step of the fit process
    (first_axis_table_dataframe, 
    organ_dimensions_table_dataframe, 
    leaf_dynamic_parameters_table_dataframe, 
    tillering_dynamic_dataframe) = fit_adel_input_data_first(plant_number, 
                                                             cohort_probabilities, 
                                                             main_stem_leaves_number_probability_distribution, 
                                                             bolting_date, 
                                                             flowering_date, 
                                                             final_axes_number) 
    
    # 3. complete the user_leaf_dynamic_parameters table
    if user_leaf_dynamic_parameters_completeness == DataCompleteness.MIN:
        leaf_dynamic_parameters_table_dataframe.ix[0]['TT_col_break'] = TT_col_break
        leaf_dynamic_parameters_table_dataframe.ix[0]['a_cohort'] = user_leaf_dynamic_parameters['a_cohort']
        leaf_dynamic_parameters_table_dataframe.ix[0]['TT_col_0'] = user_leaf_dynamic_parameters['TT_col_0']
        TT_col_nff_keys = user_leaf_dynamic_parameters['TT_col_nff'].keys()
        TT_col_nff_keys.sort()
        first_TT_col_nff = user_leaf_dynamic_parameters['TT_col_nff'][TT_col_nff_keys[0]]
        for N_cohort, leaf_dynamic_parameters_table_dataframe_grouped_by_N_cohort in leaf_dynamic_parameters_table_dataframe.groupby('N_cohort'):
            N_cohort_int = int(N_cohort)
            current_TT_col_nff = user_leaf_dynamic_parameters['TT_col_nff'][str(N_cohort_int)]
            current_dTT_MS_cohort = current_TT_col_nff - first_TT_col_nff
            index_to_set = leaf_dynamic_parameters_table_dataframe_grouped_by_N_cohort.index[0]
            leaf_dynamic_parameters_table_dataframe['TT_col_nff'][index_to_set] = current_TT_col_nff
            leaf_dynamic_parameters_table_dataframe['dTT_MS_cohort'][index_to_set] = current_dTT_MS_cohort
        leaf_dynamic_parameters_table_dataframe.ix[0]['n0'] = user_leaf_dynamic_parameters['n0']
        leaf_dynamic_parameters_table_dataframe.ix[0]['n1'] = user_leaf_dynamic_parameters['n1']
        leaf_dynamic_parameters_table_dataframe.ix[0]['n2'] = user_leaf_dynamic_parameters['n2']
    elif user_leaf_dynamic_parameters_completeness == DataCompleteness.SHORT:
        user_grouped = user_leaf_dynamic_parameters.groupby('N_cohort')
        for N_cohort, generated_group in leaf_dynamic_parameters_table_dataframe.groupby('N_cohort'):
            if N_cohort not in user_grouped.groups:
                continue
            user_group = user_grouped.get_group(N_cohort)
            index_to_get = user_group.index[0]
            index_to_set = generated_group.index[0]
            if N_cohort == 1.0:
                first_TT_col_nff = user_group['TT_col_nff'][index_to_get]
            columns_to_set = user_leaf_dynamic_parameters.columns
            leaf_dynamic_parameters_table_dataframe.ix[index_to_set][columns_to_set] = user_leaf_dynamic_parameters.ix[index_to_get]
            leaf_dynamic_parameters_table_dataframe.ix[index_to_set]['TT_col_break'] = TT_col_break
            current_TT_col_nff = user_group['TT_col_nff'][index_to_get]
            current_dTT_MS_cohort = current_TT_col_nff - first_TT_col_nff
            leaf_dynamic_parameters_table_dataframe.ix[index_to_set]['dTT_MS_cohort'] = current_dTT_MS_cohort
    elif user_leaf_dynamic_parameters_completeness == DataCompleteness.FULL:
        user_grouped_N_cohort = user_leaf_dynamic_parameters.groupby('N_cohort')
        for N_cohort, N_cohort_generated_group in leaf_dynamic_parameters_table_dataframe.groupby('N_cohort'):
            if N_cohort not in user_grouped_N_cohort.groups:
                continue
            user_N_cohort_group = user_grouped_N_cohort.get_group(N_cohort)
            N_cohort_index_to_get = user_N_cohort_group.index[0]
            N_cohort_index_to_set = N_cohort_generated_group.index[0]
            if N_cohort == 1.0:
                first_TT_col_nff = user_N_cohort_group['TT_col_nff'][N_cohort_index_to_get]
            current_TT_col_nff = user_N_cohort_group['TT_col_nff'][N_cohort_index_to_get]
            most_frequent_axis_dTT_MS_cohort = current_TT_col_nff - first_TT_col_nff
            leaf_dynamic_parameters_table_dataframe.ix[N_cohort_index_to_set]['dTT_MS_cohort'] = most_frequent_axis_dTT_MS_cohort
            highest_frequency = N_cohort_generated_group['frequency'][N_cohort_index_to_set]
            user_grouped_Nff = user_N_cohort_group.groupby('Nff')
            for Nff, Nff_generated_group in N_cohort_generated_group.groupby('Nff'):
                if Nff not in user_grouped_Nff.groups:
                    continue
                Nff_index_to_get = user_grouped_Nff.get_group(Nff).index[0]
                Nff_index_to_set = Nff_generated_group.index[0]
                if Nff_generated_group['frequency'][Nff_index_to_set] != highest_frequency:
                    current_dTT_MS_cohort = most_frequent_axis_dTT_MS_cohort + (Nff_generated_group['id_axis'][Nff_index_to_set] - N_cohort_generated_group['id_axis'][N_cohort_index_to_set]) / (4 * user_leaf_dynamic_parameters['a_cohort'][0]) 
                    leaf_dynamic_parameters_table_dataframe.ix[Nff_index_to_set]['dTT_MS_cohort'] = current_dTT_MS_cohort                   
                columns_to_set = user_leaf_dynamic_parameters.columns
                leaf_dynamic_parameters_table_dataframe.ix[Nff_index_to_set][columns_to_set] = user_leaf_dynamic_parameters.ix[Nff_index_to_get]
                leaf_dynamic_parameters_table_dataframe.ix[Nff_index_to_set]['TT_col_break'] = TT_col_break
    else:
        raise Exception('''%s is an unknown user data completeness value. \
The values can be one of %s''', (str(user_leaf_dynamic_parameters_completeness), 
                                 str([DataCompleteness.MIN, 
                                      DataCompleteness.SHORT, 
                                      DataCompleteness.FULL])))
    
    # 4. complete the user_dims table
    organ_dim_list = ['L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
    if user_dims_completeness == DataCompleteness.MIN:
        index_to_set = organ_dimensions_table_dataframe[organ_dimensions_table_dataframe['id_dim']==organ_dimensions_table_dataframe['id_dim'][0]].index
        for organ_dim in organ_dim_list:
            organ_dimensions_table_dataframe[organ_dim][index_to_set] = user_dims[organ_dim][index_to_set]
    elif user_dims_completeness == DataCompleteness.SHORT:
        user_grouped = user_dims.groupby('id_axis')
        for generated_id_dim, generated_group in organ_dimensions_table_dataframe.groupby('id_dim'):
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
                organ_dimensions_table_dataframe[organ_dim][index_to_set] = user_dims[organ_dim][index_to_get]     
    elif user_dims_completeness == DataCompleteness.FULL:
        user_grouped = user_dims.groupby('id_dim')
        for generated_id_dim, generated_group in organ_dimensions_table_dataframe.groupby('id_dim'):
            if generated_id_dim not in user_grouped.groups:
                continue
            user_group = user_grouped.get_group(generated_id_dim)
            organ_dimensions_table_dataframe.ix[generated_group.index] = user_dims.ix[user_group.index]
    else:
        raise Exception('''%s is an unknown user data completeness value. \
The values can be one of %s''', (str(user_dims_completeness), 
                                 str([DataCompleteness.MIN, 
                                      DataCompleteness.SHORT, 
                                      DataCompleteness.FULL]))) 
    
    # 5. second step of the fit process    
    (second_axis_table_dataframe, 
     absolute_second_phen_table_dataframe, 
     relative_second_phen_table_dataframe,
     absolute_organ_dimensions_table_dataframe, 
     second_leaf_dynamic_parameters_table_dataframe, 
     first_leaf_phen_table_dataframe,
     HS_GL_SSI_dynamic_dataframe, 
     relative_organ_dimensions_table_dataframe) = fit_adel_input_data_second(first_axis_table_dataframe, 
                                                                organ_dimensions_table_dataframe, 
                                                                leaf_dynamic_parameters_table_dataframe, 
                                                                GL_number, 
                                                                bolting_date, 
                                                                flowering_date, 
                                                                delais_TT_stop_del_axis, 
                                                                final_axes_number)
    
    return second_axis_table_dataframe, absolute_second_phen_table_dataframe, relative_second_phen_table_dataframe, \
           absolute_organ_dimensions_table_dataframe, second_leaf_dynamic_parameters_table_dataframe, first_leaf_phen_table_dataframe, \
           HS_GL_SSI_dynamic_dataframe, relative_organ_dimensions_table_dataframe, tillering_dynamic_dataframe
    
    

