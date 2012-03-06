'''
This module provides functions to calculate DimTable.

Created on Feb 27, 2012

@author: cchambon
'''

import pandas
import numpy


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
        
    :Types:
        - `first_parameters_table_dataframe` : pandas.DataFrame
        
    :return: The first dim table.
    :rtype: pandas.DataFrame
    '''
    id_dim_list = _create_id_dim_list(first_parameters_table_dataframe)
    absolute_index_phytomer_list = _create_absolute_index_phytomer_list(id_dim_list)
    L_blade_list = [numpy.nan for i in range(len(id_dim_list))]
    W_blade_list = [numpy.nan for i in range(len(id_dim_list))]
    L_sheath_list = [numpy.nan for i in range(len(id_dim_list))]
    W_sheath_list = [numpy.nan for i in range(len(id_dim_list))]
    L_internode_list = [numpy.nan for i in range(len(id_dim_list))]
    W_internode_list = [numpy.nan for i in range(len(id_dim_list))]
    dim_array = numpy.array([id_dim_list, absolute_index_phytomer_list, L_blade_list, W_blade_list, L_sheath_list, W_sheath_list, L_internode_list, W_internode_list]).transpose()
    return pandas.DataFrame(dim_array, columns=['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'], dtype=float)
    

def fit_dim_table_second(user_dim_table_dataframe, first_absolute_phen_table_dataframe):
    '''
    Fit the dim table: second step.
    :Parameters:
        - `user_dim_table_dataframe` : the first dim table.
        - `first_absolute_phen_table_dataframe` : the second absolute phen table.

    :Types:
        - `user_dim_table_dataframe` : int.
        - `first_absolute_phen_table_dataframe` : dict.
        
    :return: The second dim table.
    :rtype: pandas.DataFrame
    '''
    assert user_dim_table_dataframe['id_dim'].count() == user_dim_table_dataframe['id_dim'].size
    assert user_dim_table_dataframe['index_phytomer'].count() == user_dim_table_dataframe['index_phytomer'].size
    first_axis_rows_number = int(str(int(user_dim_table_dataframe['id_dim'][0]))[-2:])
    for column_name in user_dim_table_dataframe:
        assert user_dim_table_dataframe[column_name][:first_axis_rows_number].map(lambda x: x == 0.0 and 1.0 or x).fillna(0.0).all()
    dim_table_copy_dataframe = user_dim_table_dataframe.copy()
    dim_table_copy_dataframe['TT_em_phytomer'] = pandas.Series(first_absolute_phen_table_dataframe.drop(first_absolute_phen_table_dataframe.groupby('index_phytomer').groups[0.0])['TT_em_phytomer'].values)
    dim_table_copy_dataframe['L_blade'] = _fit_length(dim_table_copy_dataframe['TT_em_phytomer'], dim_table_copy_dataframe['L_blade'], first_axis_rows_number)
    dim_table_copy_dataframe['L_sheath'] = _fit_length(dim_table_copy_dataframe['TT_em_phytomer'], dim_table_copy_dataframe['L_sheath'], first_axis_rows_number)
    dim_table_copy_dataframe['L_internode'] = _fit_length(dim_table_copy_dataframe['TT_em_phytomer'], dim_table_copy_dataframe['L_internode'], first_axis_rows_number)
    
    first_axis_W_blade_tuple = (dim_table_copy_dataframe['W_blade'][0], dim_table_copy_dataframe['W_blade'][first_axis_rows_number - 1])
    first_axis_W_sheath_tuple = (dim_table_copy_dataframe['W_sheath'][0], dim_table_copy_dataframe['W_sheath'][first_axis_rows_number - 1])
    first_axis_W_internode_tuple = (dim_table_copy_dataframe['W_internode'][0], dim_table_copy_dataframe['W_internode'][first_axis_rows_number - 1])
    for name, group in dim_table_copy_dataframe[first_axis_rows_number:].groupby(by='id_dim'):
        dim_table_copy_dataframe['W_blade'][group.index] = _fit_weight(group['TT_em_phytomer'], first_axis_W_blade_tuple)
        dim_table_copy_dataframe['W_sheath'][group.index] = _fit_weight(group['TT_em_phytomer'], first_axis_W_sheath_tuple)
        dim_table_copy_dataframe['W_internode'][group.index] = _fit_weight(group['TT_em_phytomer'], first_axis_W_internode_tuple)
    
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


def _fit_length(user_dim_table_TT_em_phytomer_series, user_dim_table_L_series, first_axis_rows_number):
    '''
    Fit L_series.
    :Parameters:
        - `user_dim_table_TT_em_phytomer_series` : The TT_em_phytomer data.
        - `user_dim_table_L_series` : The lengths to fit.
        - `first_axis_rows_number` : The number of rows of the first axis. These rows are used to fit the other rows. 
    :Types:
        - `user_dim_table_TT_em_phytomer_series` : pandas.Series
        - `user_dim_table_L_series` : pandas.Series
        - `first_axis_rows_number` : int
        
    :return: The fitted lengths.
    :rtype: pandas.Series
    '''
    polynomial_coefficient_array = numpy.polyfit(user_dim_table_TT_em_phytomer_series[:first_axis_rows_number].values, user_dim_table_L_series[:first_axis_rows_number].values, 6)
    fitted_length_series = user_dim_table_L_series.copy()
    fitted_length_series[first_axis_rows_number:] = numpy.polyval(polynomial_coefficient_array, user_dim_table_TT_em_phytomer_series[first_axis_rows_number:])
    return fitted_length_series.clip_lower(0.0)


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
    x1 = user_dim_table_TT_em_phytomer_sub_series[0]
    x2 = user_dim_table_TT_em_phytomer_sub_series[-1]
    y1 = first_axis_W_tuple[0]
    y2 = first_axis_W_tuple[1]
    polynomial_coefficient_array = numpy.polyfit(numpy.array([x1, x2]), numpy.array([y1, y2]), 1)
    fitted_weight_series = pandas.Series(numpy.polyval(polynomial_coefficient_array, user_dim_table_TT_em_phytomer_sub_series))
    return fitted_weight_series.clip_lower(0.0)
    
    