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


def fit_adel_input_data(user_parameters,
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
    # 1. check input parameters validity
    assert user_parameters_completeness in DataCompleteness.__dict__.values()
    assert user_dims_completeness in DataCompleteness.__dict__.values()
    assert isinstance(TT_col_break, float)
    if user_parameters_completeness == DataCompleteness.MIN:
        # check user_parameters validity
        expected_user_parameters_keys_value_types = {'a_cohort': float, 
                                                     'TT_col_0': float, 
                                                     'TT_col_nff': float, 
                                                     'n0': float,
                                                     'n1': float,
                                                     'n2': float,
                                                     'dTT_MS_cohort': dict}
        user_parameters_keys_value_types = dict(zip(user_parameters.keys(), 
                                                    [type(value) for value in user_parameters.values()]))
        assert expected_user_parameters_keys_value_types == user_parameters_keys_value_types
        # check user_dims validity
        assert isinstance(user_dims, pandas.DataFrame)
        expected_user_dims_columns = ['index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        assert (user_dims.columns == expected_user_dims_columns).all()
        assert user_dims.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        assert user_dims['index_phytomer'].unique().size == user_dims['index_phytomer'].size
    elif user_parameters_completeness == DataCompleteness.SHORT:
        # check user_parameters validity
        assert isinstance(user_parameters, pandas.DataFrame)
        expected_user_parameters_columns = ['N_cohort', 'a_cohort', 'TT_col_0', 'TT_col_nff', 'dTT_MS_cohort', 'n0', 'n1', 'n2']
        assert (user_parameters.columns == expected_user_parameters_columns).all()
        assert user_parameters.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        assert user_parameters['N_cohort'].unique().size == user_parameters['N_cohort'].size
        # check user_dims validity
        assert isinstance(user_dims, pandas.DataFrame)
        expected_user_dims_columns = ['id_axis', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        assert (user_dims.columns == expected_user_dims_columns).all()
        assert user_dims.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        grouped = user_dims.groupby(['id_axis', 'index_phytomer'])
        assert len(grouped.groups) == user_dims.index.size            
    elif user_parameters_completeness == DataCompleteness.FULL:
        # check user_parameters validity
        assert isinstance(user_parameters, pandas.DataFrame)
        expected_user_parameters_columns = ['N_cohort', 'Nff', 'a_cohort', 'TT_col_0', 'TT_col_nff', 'dTT_MS_cohort', 'n0', 'n1', 'n2']
        assert (user_parameters.columns == expected_user_parameters_columns).all()
        assert user_parameters.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all()
        grouped = user_parameters.groupby(['N_cohort', 'Nff'])
        assert len(grouped.groups) == user_parameters.index.size    
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
    dim_table_dataframe, 
    leaf_dynamic_parameters_table_dataframe, 
    tillering_dynamic_dataframe) = fit_adel_input_data_first(plant_number, 
                                                             cohort_probabilities, 
                                                             main_stem_leaves_number_probability_distribution, 
                                                             bolting_date, 
                                                             flowering_date, 
                                                             final_axes_number) 
    
    # 3. complete the user_parameters table
    leaf_dynamic_parameters_table_dataframe.ix[0]['TT_col_break'] = TT_col_break   
    if user_parameters_completeness == DataCompleteness.MIN:
        leaf_dynamic_parameters_table_dataframe.ix[0]['a_cohort'] = user_parameters['a_cohort']
        leaf_dynamic_parameters_table_dataframe.ix[0]['TT_col_0'] = user_parameters['TT_col_0']
        leaf_dynamic_parameters_table_dataframe.ix[0]['TT_col_nff'] = user_parameters['TT_col_nff']
        for N_cohort, leaf_dynamic_parameters_table_dataframe_grouped_by_N_cohort in leaf_dynamic_parameters_table_dataframe.groupby('N_cohort'):
            N_cohort_int = int(N_cohort)
            if N_cohort_int == 1:
                current_dTT_MS_cohort = 0.0
            else:
                current_dTT_MS_cohort = user_parameters['dTT_MS_cohort'][str(N_cohort_int)]
            index_to_set = leaf_dynamic_parameters_table_dataframe_grouped_by_N_cohort.index[0]
            leaf_dynamic_parameters_table_dataframe['dTT_MS_cohort'][index_to_set] = current_dTT_MS_cohort        
        leaf_dynamic_parameters_table_dataframe.ix[0]['n0'] = user_parameters['n0']
        leaf_dynamic_parameters_table_dataframe.ix[0]['n1'] = user_parameters['n1']
        leaf_dynamic_parameters_table_dataframe.ix[0]['n2'] = user_parameters['n2']
    elif user_parameters_completeness == DataCompleteness.SHORT:
        user_grouped = user_parameters.groupby('N_cohort')
        for N_cohort, generated_group in leaf_dynamic_parameters_table_dataframe.groupby('N_cohort'):
            if N_cohort not in user_grouped.groups:
                continue
            index_to_get = user_grouped.get_group(N_cohort).index[0]
            index_to_set = generated_group.index[0]
            columns_to_set = user_parameters.columns
            leaf_dynamic_parameters_table_dataframe.ix[index_to_set][columns_to_set] = user_parameters.ix[index_to_get]        
    elif user_parameters_completeness == DataCompleteness.FULL:
        user_grouped = user_parameters.groupby(['N_cohort', 'Nff'])
        for (N_cohort, Nff), generated_group in leaf_dynamic_parameters_table_dataframe.groupby(['N_cohort', 'Nff']):
            if (N_cohort, Nff) not in user_grouped.groups:
                continue
            index_to_get = user_grouped.get_group((N_cohort, Nff)).index[0]
            index_to_set = generated_group.index[0]
            columns_to_set = user_parameters.columns
            leaf_dynamic_parameters_table_dataframe.ix[index_to_set][columns_to_set] = user_parameters.ix[index_to_get]
    else:
        raise Exception('''%s is an unknown user data completeness value. \
The values can be one of %s''', (str(user_parameters_completeness), 
                                 str([DataCompleteness.MIN, 
                                      DataCompleteness.SHORT, 
                                      DataCompleteness.FULL])))
    
    # 4. complete the user_dims table
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
    
    # 5. second step of the fit process    
    (second_axis_table_dataframe, 
     absolute_second_phen_table_dataframe, 
     relative_second_phen_table_dataframe,
     absolute_dim_table_dataframe, 
     second_leaf_dynamic_parameters_table_dataframe, 
     first_leaf_phen_table_dataframe,
     HS_GL_SSI_dynamic_dataframe, 
     relative_dim_table_dataframe) = fit_adel_input_data_second(first_axis_table_dataframe, 
                                                                dim_table_dataframe, 
                                                                leaf_dynamic_parameters_table_dataframe, 
                                                                GL_number, 
                                                                bolting_date, 
                                                                flowering_date, 
                                                                delais_TT_stop_del_axis, 
                                                                final_axes_number)
    
    return second_axis_table_dataframe, absolute_second_phen_table_dataframe, relative_second_phen_table_dataframe, \
           absolute_dim_table_dataframe, second_leaf_dynamic_parameters_table_dataframe, first_leaf_phen_table_dataframe, \
           HS_GL_SSI_dynamic_dataframe, relative_dim_table_dataframe, tillering_dynamic_dataframe
    
    

