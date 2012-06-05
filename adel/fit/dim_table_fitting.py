'''
This module provides functions to calculate DimTable.

Created on Feb 27, 2012

@author: cchambon
'''

import pandas
import numpy as np


def fit_dim_table_first(first_parameters_table_dataframe):
    '''
    Fit the dim table: first step.
    :Parameters:
        - `first_parameters_table_dataframe` : the first parameters table.
            * N_cohort: the cohort number.
            * id_axis: the identifier made from concatenation of N_cohort and Nff.
            * frequency: the occurrence frequency of id_axis.
            * Nff: the final leaves number of the current axis.
            * a_cohort: The slope of phytomer emergence.
            * TT_col_0: The Thermal Time for Haun Stage equal to 0.
            * TT_col_break: The Thermal Time when the slope of phytomers emergence is changing.
            * TT_col_nff: The Thermal Time when Haun Stage is equal to the total number of vegetative phytomers formed on an axis.
            * dTT_MS_cohort: The thermal time delta between the main stem flowering and the current raw tiller flowering.
            * n0: ???
            * n1: ???
            * n2: ???
            * t0: ???
            * t1: ???
            * t2: ???
            * hs_t1: ???
            * a: ???
            * c: ???
            * d: ???
            
    :Types:
        - `first_parameters_table_dataframe` : pandas.DataFrame
        
    :return: The first dim table.
    :rtype: pandas.DataFrame
    '''
    id_dim_list = _create_id_dim_list(first_parameters_table_dataframe)
    absolute_index_phytomer_list = _create_absolute_index_phytomer_list(id_dim_list)
    L_blade_list = [np.nan for i in range(len(id_dim_list))]
    W_blade_list = [np.nan for i in range(len(id_dim_list))]
    L_sheath_list = [np.nan for i in range(len(id_dim_list))]
    W_sheath_list = [np.nan for i in range(len(id_dim_list))]
    L_internode_list = [np.nan for i in range(len(id_dim_list))]
    W_internode_list = [np.nan for i in range(len(id_dim_list))]
    dim_array = np.array([id_dim_list, absolute_index_phytomer_list, L_blade_list, W_blade_list, L_sheath_list, W_sheath_list, L_internode_list, W_internode_list]).transpose()
    return pandas.DataFrame(dim_array, columns=['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'], dtype=float)
    

def fit_dim_table_second(user_dim_table_dataframe, absolute_second_phen_table_dataframe):
    '''
    Fit the dim table: second step.
    :Parameters:
        - `user_dim_table_dataframe` : the first dim table.
        - `absolute_second_phen_table_dataframe` : the second absolute phen table.

    :Types:
        - `user_dim_table_dataframe` : int.
        - `absolute_second_phen_table_dataframe` : dict.
        
    :return: The second dim table.
    :rtype: pandas.DataFrame
    '''
    assert user_dim_table_dataframe['id_dim'].count() == user_dim_table_dataframe['id_dim'].size
    assert user_dim_table_dataframe['index_phytomer'].count() == user_dim_table_dataframe['index_phytomer'].size
    first_axis_rows_number = int(str(int(user_dim_table_dataframe['id_dim'][0]))[-2:])
    for column_name in user_dim_table_dataframe:
        assert user_dim_table_dataframe[column_name][:first_axis_rows_number].map(lambda x: x == 0.0 and 1.0 or x).fillna(0.0).all()
    first_row_number_to_fit = user_dim_table_dataframe['L_blade'].last_valid_index() + 1
    dim_table_copy_dataframe = user_dim_table_dataframe.copy()
    dim_table_copy_dataframe['TT_em_phytomer'] = pandas.Series(absolute_second_phen_table_dataframe.drop(absolute_second_phen_table_dataframe.groupby('index_phytomer').groups[0.0])['TT_em_phytomer'].values)
    dim_table_copy_dataframe['L_blade'] = _fit_length(dim_table_copy_dataframe['TT_em_phytomer'], dim_table_copy_dataframe['L_blade'], first_axis_rows_number, first_row_number_to_fit)
    dim_table_copy_dataframe['L_sheath'] = _fit_length(dim_table_copy_dataframe['TT_em_phytomer'], dim_table_copy_dataframe['L_sheath'], first_axis_rows_number, first_row_number_to_fit)
    dim_table_copy_dataframe['L_internode'] = _fit_length(dim_table_copy_dataframe['TT_em_phytomer'], dim_table_copy_dataframe['L_internode'], first_axis_rows_number, first_row_number_to_fit)
    
    first_axis_W_blade_tuple = (dim_table_copy_dataframe['W_blade'][0], dim_table_copy_dataframe['W_blade'][first_axis_rows_number - 1])
    first_axis_W_sheath_tuple = (dim_table_copy_dataframe['W_sheath'][0], dim_table_copy_dataframe['W_sheath'][first_axis_rows_number - 1])
    first_axis_W_internode_tuple = (dim_table_copy_dataframe['W_internode'][0], dim_table_copy_dataframe['W_internode'][first_axis_rows_number - 1])
    for name, group in dim_table_copy_dataframe[first_row_number_to_fit:].groupby(by='id_dim'):
        dim_table_copy_dataframe['W_blade'][group.index] = _fit_weight(group['TT_em_phytomer'], first_axis_W_blade_tuple)
        dim_table_copy_dataframe['W_sheath'][group.index] = _fit_weight(group['TT_em_phytomer'], first_axis_W_sheath_tuple)
        dim_table_copy_dataframe['W_internode'][group.index] = _fit_weight(group['TT_em_phytomer'], first_axis_W_internode_tuple)
    
#    # uncomment for dimTable validation
#    for column_name in ['L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']:
#        from matplotlib import pyplot
#        from openalea.core.path import path
#        pyplot.figure()
#        legend_list = []
#        for id_dim, group in dim_table_copy_dataframe.groupby(by='id_dim'):
#            current_column_array = group[column_name].values
#            legend_list.append(str(int(id_dim)))
#            pyplot.plot(current_column_array, group['index_phytomer'].values)
#
#        pyplot.xlabel(column_name)
#        pyplot.ylabel('index_phytomer')
#        pyplot.title(column_name)
#        pyplot.legend(legend_list, 'best')
#        pyplot.savefig(path('data/test_fitting2/dimTable_validation/%s.png') % column_name)
    
    return dim_table_copy_dataframe.drop(['TT_em_phytomer'], axis=1)


def _create_id_dim_list(first_parameters_dataframe):
    '''
    Create id_dim list.
    :Parameters:
        - `first_parameters_dataframe` : the first observations table, with the following headers: 
            * N_cohort: the cohort number.
            * id_axis: the identifier made from concatenation of N_cohort and Nff.
            * frequency: the occurrence frequency of id_axis.
            * Nff: the final leaves number of the current axis.
            * a_cohort: The slope of phytomer emergence.
            * TT_col_0: The Thermal Time for Haun Stage equal to 0.
            * TT_col_break: The Thermal Time when the slope of phytomers emergence is changing.
            * TT_col_nff: The Thermal Time when Haun Stage is equal to the total number of vegetative phytomers formed on an axis.
            * dTT_MS_cohort: The thermal time delta between the main stem flowering and the current raw tiller flowering.
            * n0: ???
            * n1: ???
            * n2: ???
            * t0: ???
            * t1: ???
            * t2: ???
            * hs_t1: ???
            * a: ???
            * c: ???
            * d: ???
            
    :Types:
        - `first_parameter_dataframe` : pandas.DataFrame
        
    :return: The new completed id_dim list.
    :rtype: list
    '''
    sorted_id_axis = first_parameters_dataframe['id_axis']
    id_dim_list = []
    for id_dim in sorted_id_axis:
        N_phyt = int(str(int(id_dim))[-2:])
        for i in range(N_phyt):
            id_dim_list.append(id_dim)
    return id_dim_list


def _create_absolute_index_phytomer_list(first_parameters_table_id_dim_list):
    '''
    Create list of absolute phytomer index.
    :Parameters:
        - `id_dim_list` : the id_dim list.
    :Types:
        - `id_dim_list` : list
        
    :return: The absolute_index_phytomer list.
    :rtype: list
    '''
    absolute_index_phytomer_list = []
    i = 0
    while i < len(first_parameters_table_id_dim_list):
        N_phyt = int(str(int(first_parameters_table_id_dim_list[i]))[-2:])
        for j in range(1, N_phyt + 1):
            absolute_index_phytomer_list.append(j)    
        i = i + j
    return absolute_index_phytomer_list


def _fit_length(user_dim_table_TT_em_phytomer_series, user_dim_table_L_series, first_axis_rows_number, first_row_number_to_fit):
    '''
    Fit L_series.
    :Parameters:
        - `user_dim_table_TT_em_phytomer_series` : The TT_em_phytomer data.
        - `user_dim_table_L_series` : The lengths to fit.
        - `first_axis_rows_number` : The number of rows of the first axis. These rows are used to fit the other rows. 
        - `first_row_number_to_fit` : The number of the first row to fit. 
    :Types:
        - `user_dim_table_TT_em_phytomer_series` : pandas.Series
        - `user_dim_table_L_series` : pandas.Series
        - `first_axis_rows_number` : int
        - `first_row_number_to_fit` : int 
        
    :return: The fitted lengths.
    :rtype: pandas.Series
    '''
    polynomial_coefficient_array = np.polyfit(user_dim_table_TT_em_phytomer_series[:first_axis_rows_number].values, user_dim_table_L_series[:first_axis_rows_number].values, 6)
    fitted_length_series = user_dim_table_L_series.copy()
    fitted_length_series[first_row_number_to_fit:] = np.polyval(polynomial_coefficient_array, user_dim_table_TT_em_phytomer_series[first_row_number_to_fit:])
    main_stem_last_TT_em_phytomer = user_dim_table_TT_em_phytomer_series[first_axis_rows_number - 1]
    indexes_to_change = user_dim_table_TT_em_phytomer_series[user_dim_table_TT_em_phytomer_series > main_stem_last_TT_em_phytomer].index
    indexes_to_change = indexes_to_change[indexes_to_change >= first_row_number_to_fit]
    fitted_length_series[indexes_to_change] = user_dim_table_L_series[first_axis_rows_number - 1]
    return fitted_length_series


def _fit_weight(user_dim_table_TT_em_phytomer_sub_series, first_axis_W_tuple):
    '''
    Fit weight.
    :Parameters:
        - `user_dim_table_TT_em_phytomer_sub_series` : The TT_em_phytomer data, for one axis.
        - `first_axis_W_tuple` : The first phytomer weights of the first axis and the last phytomer weights of the first axis. 
    :Types:
        - `user_dim_table_TT_em_phytomer_sub_series` : pandas.Series
        - `first_axis_W_tuple` : tuple
        
    :return: The fitted weigths.
    :rtype: pandas.Series
    '''
    x1 = user_dim_table_TT_em_phytomer_sub_series[user_dim_table_TT_em_phytomer_sub_series.index[0]]
    x2 = user_dim_table_TT_em_phytomer_sub_series[user_dim_table_TT_em_phytomer_sub_series.index[-1]]
    y1 = first_axis_W_tuple[0]
    y2 = first_axis_W_tuple[1]
    polynomial_coefficient_array = np.polyfit(np.array([x1, x2]), np.array([y1, y2]), 1)
    fitted_weight_series = pandas.Series(np.polyval(polynomial_coefficient_array, user_dim_table_TT_em_phytomer_sub_series))
    return fitted_weight_series.clip_lower(0.0)


def create_dim_table_relative_dataframe(absolute_dim_table_dataframe):
    '''
    Create the relative dim table table dataframe from the absolute one. 
    :Parameters:
        - `absolute_dim_table_dataframe` : the absolute dim table dataframe.
    :Types:
        - `absolute_dim_table_dataframe` : pandas.Dataframe
        
    :return: the relative phen table dataframe.
    :rtype: pandas.Dataframe
    '''
    assert absolute_dim_table_dataframe.count().max() == absolute_dim_table_dataframe.count().min() == absolute_dim_table_dataframe.index.size
    relative_dim_table_dataframe = pandas.DataFrame(index=absolute_dim_table_dataframe.index, columns=['id_dim', 'index_rel_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'])
    relative_dim_table_dataframe[['id_dim', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']] = absolute_dim_table_dataframe[['id_dim', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']]
    tmp_series = pandas.Series(relative_dim_table_dataframe.index)
    for name, group in absolute_dim_table_dataframe.groupby('id_dim'):
        def normalize_index_phytomer(i):
            return group['index_phytomer'][i] / float(group['index_phytomer'][group.index[-1]])
        relative_dim_table_dataframe['index_rel_phytomer'][group.index] = tmp_series[group.index].map(normalize_index_phytomer) 
        
    return relative_dim_table_dataframe
        
        
    