# -*- coding: latin-1 -*-
'''
This module provides functions to calculate PhenTable values.

Created on 28 nov. 2011

@author: cchambon
'''

import pandas

em_col_phytomer_delay = 1.6 #TODO: to change!

def create_id_phen_list(fitted_parameter_dataframe):
    '''
    Create list of id_phen.
    :Parameters:
        - `fitted_parameter_dataframe` : the fitted observations table, with the following headers: 
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
        - `fitted_parameter_dataframe` : pandas.DataFrame
        
    :return: The new completed id_phen list.
    :rtype: list
    '''
    sorted_id_axis = fitted_parameter_dataframe['id_axis']
    id_phen_list = []
    for id_phen in sorted_id_axis:
        N_phyt = int(str(int(id_phen))[-2:])
        for i in range(N_phyt + 1):
            id_phen_list.append(id_phen)
    return id_phen_list


def create_absolute_index_phytomer_list(id_phen_list):
    '''
    Create list of absolute phytomer index.
    :Parameters:
        - `id_phen_list` : the id_phen list.
    :Types:
        - `id_phen_list` : list
        
    :return: The absolute_index_phytomer list.
    :rtype: list
    '''
    absolute_index_phytomer_list = []
    i = 0
    while i < len(id_phen_list):
        N_phyt = int(str(int(id_phen_list[i]))[-2:])
        for j in range(N_phyt + 1):
            absolute_index_phytomer_list.append(j)    
        i = i + j + 1
    return absolute_index_phytomer_list


def create_absolute_TT_col_phytomer_list(fitted_parameter_dataframe):
    '''
    Create list of TT_col_phytomer. 
    :Parameters:
        - `fitted_parameter_dataframe` : the fitted observations table, with the following headers: 
            * N_cohort: the cohort number.
            * id_axis: the identifier made from concatenation of N_cohort and Nff.
            * frequency: the occurrence frequency of id_axis.
            * Nff: the final leaves number of the current axis.
            * a_cohort: The slope of phytomer emergence.
            * TT_col_0: The Thermal Time for Haun Stage equal to 0.
            * TT_col_break: The Thermal Time when the slope of phytomers emergence is changing.
            * TT_col_nff: The Thermal Time when Haun Stage is equal to the total number of vegetative phytomers formed on an axis.
            * dTT_MS_cohort: The thermal time delta between the main stem flowering and the current raw tiller flowering.
          The table is completely filled and ordered by frequency. 
    :Types:
        - `fitted_parameter_dataframe` : pandas.DataFrame
        
    :return: The TT_col_phytomer list.
    :rtype: list
    '''
    assert fitted_parameter_dataframe.count().max() == fitted_parameter_dataframe.count().min() == fitted_parameter_dataframe.index.size
    TT_col_phytomer_list = []
    
    if fitted_parameter_dataframe['TT_col_break'][0] == 0.0:
        for i in fitted_parameter_dataframe.index:
            TT_col_0_i = fitted_parameter_dataframe['TT_col_0'][i]
            a_cohort_i = fitted_parameter_dataframe['a_cohort'][i]
            Nff_i = fitted_parameter_dataframe['Nff'][i]
            phytomer_indexes = range(int(Nff_i) + 1)
            for j in phytomer_indexes:
                TT_col_phytomer_i_j = TT_col_0_i + j / a_cohort_i
                TT_col_phytomer_list.append(TT_col_phytomer_i_j)
    else: # bilinear calculation method
        for i in fitted_parameter_dataframe.index:
            TT_col_nff_i = fitted_parameter_dataframe['TT_col_nff'][i]
            TT_col_break_i = fitted_parameter_dataframe['TT_col_break'][i]
            Nff_i = fitted_parameter_dataframe['Nff'][i]
            a_cohort_i = fitted_parameter_dataframe['a_cohort'][i]
            TT_col_0_i = fitted_parameter_dataframe['TT_col_0'][i]
            phytomer_indexes = range(Nff_i + 1)
            HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
            a2_i = (Nff_i - HS_break_i) / (TT_col_nff_i - TT_col_break_i)
            for j in phytomer_indexes: 
                if j < HS_break_i:
                    TT_col_phytomer_i_j = TT_col_0_i + j / a_cohort_i
                else:
                    TT_col_phytomer_i_j = (j - HS_break_i) / a2_i + TT_col_break_i
                TT_col_phytomer_list.append(TT_col_phytomer_i_j) 
    
    return TT_col_phytomer_list
    
    
def create_absolute_TT_em_phytomer_list(TT_col_phytomer_list, fitted_parameter_dataframe):
    '''
    Create list of TT_em_phytomer. 
    :Parameters:
        - `TT_col_phytomer_list` : The TT_col_phytomer list.
        - `fitted_parameter_dataframe` : the fitted observations table, with the following headers: 
            * N_cohort: the cohort number.
            * id_axis: the identifier made from concatenation of N_cohort and Nff.
            * frequency: the occurrence frequency of id_axis.
            * Nff: the final leaves number of the current axis.
            * a_cohort: The slope of phytomer emergence.
            * TT_col_0: The Thermal Time for Haun Stage equal to 0.
            * TT_col_break: The Thermal Time when the slope of phytomers emergence is changing.
            * TT_col_nff: The Thermal Time when Haun Stage is equal to the total number of vegetative phytomers formed on an axis.
            * dTT_MS_cohort: The thermal time delta between the main stem flowering and the current raw tiller flowering.
          The table is completely filled and ordered by frequency. 
    :Types:
        - `TT_col_phytomer_list` : list
        - `fitted_parameter_dataframe` : pandas.DataFrame
        
    :return: The TT_em_phytomer list.
    :rtype: list
    '''
    assert fitted_parameter_dataframe.count().max() == fitted_parameter_dataframe.count().min() == fitted_parameter_dataframe.index.size
    TT_em_phytomer_list = []
    for i in fitted_parameter_dataframe.index:
        Nff_i = fitted_parameter_dataframe['Nff'][i]
        phytomer_indexes = range(int(Nff_i) + 1)
        TT_col_break_i = fitted_parameter_dataframe['TT_col_break'][i]
        a_cohort_i = fitted_parameter_dataframe['a_cohort'][i]
        for j in phytomer_indexes: 
            TT_col_phytomer_j = TT_col_phytomer_list[j]
            if TT_col_phytomer_j < TT_col_break_i:
                TT_em_phytomer_j = TT_col_phytomer_j - (em_col_phytomer_delay / a_cohort_i)
            else:
                TT_col_0_i = fitted_parameter_dataframe['TT_col_0'][i]
                HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
                TT_col_nff_i = fitted_parameter_dataframe['TT_col_nff'][i]
                a2_i = (Nff_i - HS_break_i) / (TT_col_nff_i - TT_col_break_i)
                tmp_res = TT_col_phytomer_j - (em_col_phytomer_delay / a2_i)
                if tmp_res > HS_break_i:
                    TT_em_phytomer_j = tmp_res
                else:
                    TT_em_phytomer_j = TT_col_break_i - (em_col_phytomer_delay - a2_i * (TT_col_phytomer_j - TT_col_break_i)) / a_cohort_i
            TT_em_phytomer_list.append(TT_em_phytomer_j)
            
    return TT_em_phytomer_list


def create_phen_table_relative_dataframe(absolute_phen_table_dataframe, first_leaf_phen_table_dataframe):
    '''
    Create the relative phen table table dataframe from the absolute one. 
    :Parameters:
        - `absolute_phen_table_dataframe` : the absolute phen table dataframe.
        - `first_leaf_phen_table_dataframe` : the first leaf phen table dataframe.
    :Types:
        - `absolute_phen_table_dataframe` : pandas.Dataframe
        - `first_leaf_phen_table_dataframe` : pandas.Dataframe
        
    :return: the relative phen table table dataframe.
    :rtype: pandas.Dataframe
    '''
    relative_phen_table_dataframe = absolute_phen_table_dataframe.copy()
    tmp_series = pandas.Series(relative_phen_table_dataframe.index)
    for name, group in absolute_phen_table_dataframe.groupby('id_phen'):
        def normalize_index_phytomer(i):
            return group['index_phytomer'][i] / group['index_phytomer'][-1]
        relative_phen_table_dataframe['index_phytomer'][group.index] = tmp_series[group.index].map(normalize_index_phytomer)
        def get_relative_TT_col_phytomer(i):
            return group['TT_col_phytomer'][i] - first_leaf_phen_table_dataframe['TT_col_phytomer'][name]
        relative_phen_table_dataframe['TT_col_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_col_phytomer)
        def get_relative_TT_em_phytomer(i):
            return group['TT_em_phytomer'][i] - first_leaf_phen_table_dataframe['TT_em_phytomer'][name]
        relative_phen_table_dataframe['TT_em_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_em_phytomer)
        
    return relative_phen_table_dataframe

