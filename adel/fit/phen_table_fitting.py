# -*- coding: latin-1 -*-
'''
This module provides functions to calculate PhenTable.

Created on 28 nov. 2011

@author: cchambon
'''

import numpy as np
import pandas

delais_phyll_col_tip = 1.6
delais_phyll_sen_disp = 3.0
TT_del_Fhaut = 3000
Nbr_Fhaut_persistant = 5


def fit_phen_table_second(second_parameters_table_dataframe):
    '''
    Fit the phen table: first step.
    :Parameters:
        - `second_parameters_table_dataframe` : the fitted observations table, with the following headers: 
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
            * d: ??
    :Types:
        - `second_parameters_table_dataframe` : pandas.DataFrame
        
    :return: The phen table. 
    :rtype: pandas.DataFrame
    '''
    assert second_parameters_table_dataframe.count().max() == second_parameters_table_dataframe.count().min() == second_parameters_table_dataframe.index.size
    id_phen_list = _create_id_phen_list(second_parameters_table_dataframe)
    absolute_index_phytomer_list = _create_absolute_index_phytomer_list(id_phen_list)
    absolute_TT_col_phytomer_list = _create_absolute_TT_col_phytomer_list(second_parameters_table_dataframe)
    absolute_TT_em_phytomer_list = _create_absolute_TT_em_phytomer_list(absolute_TT_col_phytomer_list, second_parameters_table_dataframe)
    absolute_TT_sen_phytomer_list = _create_absolute_TT_sen_phytomer_list(second_parameters_table_dataframe)
    absolute_TT_del_phytomer_list = _create_absolute_TT_del_phytomer_list(id_phen_list, absolute_TT_sen_phytomer_list, second_parameters_table_dataframe)
    absolute_phen_table_array = np.array([id_phen_list, absolute_index_phytomer_list, absolute_TT_em_phytomer_list, absolute_TT_col_phytomer_list, absolute_TT_sen_phytomer_list, absolute_TT_del_phytomer_list]).transpose()
    return pandas.DataFrame(absolute_phen_table_array, columns=['id_phen', 'index_phytomer', 'TT_em_phytomer', 'TT_col_phytomer', 'TT_sen_phytomer', 'TT_del_phytomer'], dtype=float)


def _create_id_phen_list(second_parameters_table_dataframe):
    '''
    Create list of id_phen.
    :Parameters:
        - `second_parameters_table_dataframe` : the fitted observations table, with the following headers: 
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
            * d: ??
    :Types:
        - `second_parameters_table_dataframe` : pandas.DataFrame
        
    :return: The new completed id_phen list.
    :rtype: list
    '''
    sorted_id_axis = second_parameters_table_dataframe['id_axis']
    id_phen_list = []
    for id_phen in sorted_id_axis:
        N_phyt = int(str(int(id_phen))[-2:])
        for i in range(N_phyt + 1):
            id_phen_list.append(id_phen)
    return id_phen_list


def _create_absolute_index_phytomer_list(second_parameters_table_id_phen_list):
    '''
    Create list of absolute phytomer index.
    :Parameters:
        - `second_parameters_table_id_phen_list` : the id_phen list.
    :Types:
        - `second_parameters_table_id_phen_list` : list
        
    :return: The absolute_index_phytomer list.
    :rtype: list
    '''
    absolute_index_phytomer_list = []
    i = 0
    while i < len(second_parameters_table_id_phen_list):
        N_phyt = int(str(int(second_parameters_table_id_phen_list[i]))[-2:])
        for j in range(N_phyt + 1):
            absolute_index_phytomer_list.append(j)    
        i = i + j + 1
    return absolute_index_phytomer_list


def _create_absolute_TT_col_phytomer_list(second_parameters_table_dataframe):
    '''
    Create list of TT_col_phytomer. 
    :Parameters:
        - `second_parameters_table_dataframe` : the fitted observations table, with the following headers: 
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
            * d: ??
          The table is completely filled and ordered by frequency. 
    :Types:
        - `second_parameters_table_dataframe` : pandas.DataFrame
        
    :return: The TT_col_phytomer list.
    :rtype: list
    '''
    TT_col_phytomer_list = []
    for i in second_parameters_table_dataframe.index:
        N_cohort_i, id_axis_i, frequency_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, dTT_MS_cohort_i, n0_i, n1_i, n2_i, t0_i, t1_i, t2_i, hs_t1_i, a_i, c_i, d_i, RMSE_gl = second_parameters_table_dataframe.ix[i].tolist()
        phytomer_indexes = range(int(Nff_i + 1))
        HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
        a2_i = (Nff_i - HS_break_i) / (TT_col_nff_i - TT_col_break_i)
        for j in phytomer_indexes: 
            if j < HS_break_i: # linear mode
                TT_col_phytomer_i_j = TT_col_0_i + j / a_cohort_i
            else: # bilinear mode
                TT_col_phytomer_i_j = (j - HS_break_i) / a2_i + TT_col_break_i
            TT_col_phytomer_list.append(TT_col_phytomer_i_j) 
    
    return TT_col_phytomer_list
    
    
def _create_absolute_TT_em_phytomer_list(phen_table_TT_col_phytomer_list, second_parameters_table_dataframe):
    '''
    Create list of TT_em_phytomer. 
    :Parameters:
        - `phen_table_TT_col_phytomer_list` : The TT_col_phytomer list.
        - `second_parameters_table_dataframe` : the fitted observations table, with the following headers: 
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
            * d: ??
          The table is completely filled and ordered by frequency. 
    :Types:
        - `phen_table_TT_col_phytomer_list` : list
        - `second_parameters_table_dataframe` : pandas.DataFrame
        
    :return: The TT_em_phytomer list.
    :rtype: list
    '''
    TT_em_phytomer_list = []
    current_TT_em_phytomer_row_index = 0
    for i in second_parameters_table_dataframe.index:
        N_cohort_i, id_axis_i, frequency_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, dTT_MS_cohort_i, n0_i, n1_i, n2_i, t0_i, t1_i, t2_i, hs_t1_i, a_i, c_i, d_i, RMSE_gl = second_parameters_table_dataframe.ix[i].tolist()
        phytomer_indexes = np.arange(int(Nff_i) + 1) + current_TT_em_phytomer_row_index
        j = 0
        for j in phytomer_indexes: 
            TT_col_phytomer_j = phen_table_TT_col_phytomer_list[j]
            if TT_col_phytomer_j < TT_col_break_i:
                TT_em_phytomer_j = TT_col_phytomer_j - (delais_phyll_col_tip / a_cohort_i)
            else:
                HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
                a2_i = (Nff_i - HS_break_i) / (TT_col_nff_i - TT_col_break_i)
                tmp_res = TT_col_phytomer_j - (delais_phyll_col_tip / a2_i)
                if tmp_res > HS_break_i:
                    TT_em_phytomer_j = tmp_res
                else:
                    TT_em_phytomer_j = TT_col_break_i - (delais_phyll_col_tip - a2_i * (TT_col_phytomer_j - TT_col_break_i)) / a_cohort_i
            TT_em_phytomer_list.append(TT_em_phytomer_j)
        current_TT_em_phytomer_row_index += phytomer_indexes.size
            
    return TT_em_phytomer_list


def _create_absolute_TT_sen_phytomer_list(second_parameters_table_dataframe):
    '''
    Create list of TT_sen_phytomer. To compute TT_sen_phytomer values, the algorithm main steps are:
        - calculate HS(TT):
            HS_break[i] = a_cohort[i] * (TT_col_break[i] - TT_col_0[i])
            a2[i] = (Nff[i] - HS_break[i]) / (TT_col_nff[i] - TT_col_break[i])
            if j < HS_break[i]: # linear mode
                if a_cohort[i]*(TT-TT_col_0[i])<Nff[i] then HS(TT) = a_cohort[i]*(TT-TT_col_0[i])
                else HS(TT) = Nff[i]
            else: # bilinear mode
                if (TT - TT_col_break[i]) * a2[i] + HS_break[i]<Nff[i] then HS(TT) = (TT - TT_col_break[i]) * a2[i] + HS_break[i]
                else HS(TT) = Nff[i]
        - calculate GL(TT):
            - for the main stem:
              if TT<t0[i] then GL(TT) = HS(TT)
              else if TT<t1[i] then GL(TT) = n0[i]+(n1[i]-n0[i])*(TT-t0[i])/(t1[i]-t0[i])
                   else if TT<t2[i] then GL(TT) = n1[i]+(n2[i]-n1[i])*(TT-t1[i])/(t2[i]-t1[i])
                        else if a[i]*(TT-t2[i])^3+c[i]*(TT-t2[i])+d[i]>0 then GL(TT) = a[i]*(TT-t2[i])^3+c[i]*(TT-t2[i])+d[i]
                             else GL(TT) = 0
            - for the tillers:
              if TT<t0[i] then GL(TT) = HS(TT)
              else if TT<t1[i] then GL(TT) = n0[i]
                   else if TT<=t2[i] then GL(TT) = n0[i]+(n2[i]-n0[i])*(TT-t1[i])/(t2[i]-t1[i])
                        else if a[i]*(TT-t2[i])^3+c[i]*(TT-t2[i])+d[i]>0 then GL(TT) = a[i]*(TT-t2[i])^3+c[i]*(TT-t2[i])+d[i]
                             else GL(TT) = 0
        - calculate SSI(TT):
            SSI(TT) = HS(TT) - GL(TT)
        - find TT_sen_phytomer so that j = SSI(TT_sen_phytomer), i.e. SSI(TT_sen_phytomer) - j = 0
        In all this algorithm, TT is the linear thermal time, i is an axis identifier and j is a phytomer of i.
    :Parameters:
        - `second_parameters_table_dataframe` : the fitted observations table, with the following headers:
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
          The table is completely filled and ordered by frequency. 
    :Types:
        - `second_parameters_table_dataframe` : pandas.DataFrame
        
    :return: The TT_em_phytomer list.
    :rtype: list
    '''
    TT_sen_phytomer_list = []
    
    def get_real_roots(poly):
        roots_array = poly.r
        null_imag_solutions = filter(lambda x: x.imag == 0.0, roots_array)
        real_solutions = map(lambda x: x.real, null_imag_solutions)
        return np.array(real_solutions)
        
    for i in second_parameters_table_dataframe.index:
        TT_sen_phytomer_i_list = []        
        N_cohort_i, id_axis_i, frequency_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, dTT_MS_cohort_i, n0_i, n1_i, n2_i, t0_i, t1_i, t2_i, hs_t1_i, a_i, c_i, d_i, RMSE_gl = second_parameters_table_dataframe.ix[i].tolist()
        HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
        HS_1, HS_2, GL_2, GL_3, GL_4 = _create_HS_GL_polynomial(HS_break_i, N_cohort_i, id_axis_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, n0_i, n1_i, n2_i, t0_i, t1_i, t2_i, a_i, c_i, d_i)
        phytomer_indexes_i = np.arange(int(Nff_i) + 1)
#        # Uncomment to create the plots from roots calculation
#        TT_to_plot_method_2 = []
#        HS_to_plot_method_2 = []
#        GL_to_plot_method_2 = []
#        SSI_to_plot_method_2 = []

        for j in phytomer_indexes_i:
            # define HS according to j
            if j < HS_break_i: # linear mode
                HS = HS_1
            else: # bilinear mode
                HS = HS_2
            GL_1 = HS
            # Suppose we are in the [0,t0_i] phase.
            GL = GL_1
            SSI = HS - GL
            if j == 0:
                TT_sen_phytomer_i_list.append(t0_i)
            else:
                # Find (SSI - j) real root.
                SSI_root_array = get_real_roots(SSI - j)
                if SSI_root_array.size == 0 or SSI_root_array[0] > t0_i:
                    # Suppose we are in the ]t0_i,t1_i] phase.
                    GL = GL_2
                    SSI = HS - GL
                    # Find (SSI - j) real root.
                    SSI_root_array = get_real_roots(SSI - j)
                    if SSI_root_array.size == 0 or SSI_root_array[0] <= t0_i or SSI_root_array[0] > t1_i:
                        # Suppose we are in the ]t1_i,t2_i] phase.
                        GL = GL_3
                        SSI = HS - GL
                        # Find (SSI - j) real root.
                        SSI_root_array = get_real_roots(SSI - j)
                        if SSI_root_array.size == 0 or SSI_root_array[0] <= t1_i or SSI_root_array[0] > t2_i:
                            # We must be in the ]t2_i,infinity[ phase.
                            GL = GL_4
                            SSI = HS - GL
                            # Find (SSI - j) real root.
                            SSI_root_array = get_real_roots(SSI - j)
                            if SSI_root_array.size == 0 or SSI_root_array[0] <= t2_i:
                                raise Exception('ERROR !!!!! This shouldn\'t occurred')
                            if HS(SSI_root_array[0]) > Nff_i:
                                HS = np.poly1d([Nff_i])
                            if GL(SSI_root_array[0]) < 0.0:
                                GL = np.poly1d([0.0])
                            SSI = HS - GL
                            # Find (SSI - j) real root again.
                            SSI_root_array = get_real_roots(SSI - j)
                            if SSI_root_array.size == 0 or SSI_root_array[0] <= t2_i:
                                raise Exception('ERROR !!!!! This shouldn\'t occurred')   
                TT_sen_phytomer_i_list.append(SSI_root_array[0])
                
#            # Uncomment to create the plots from roots calculation
#            TT_to_plot_method_2.append(TT_sen_phytomer_i_list[-1])
#            HS_to_plot_method_2.append(HS(TT_sen_phytomer_i_list[-1]))
#            GL_to_plot_method_2.append(GL(TT_sen_phytomer_i_list[-1]))
#            SSI_to_plot_method_2.append(j)
#        # Uncomment to create the plots
#        from matplotlib import pyplot
#        from openalea.core.path import path
#        pyplot.figure(int(id_axis_i))
#        pyplot.plot(TT_to_plot_method_2, HS_to_plot_method_2, 'b', marker='o')
#        pyplot.plot(TT_to_plot_method_2, GL_to_plot_method_2, 'g', marker='o')
#        pyplot.plot(TT_to_plot_method_2, SSI_to_plot_method_2, 'r', marker='o')
#        pyplot.xlabel('TT')
#        pyplot.ylabel('Phytomer')
#        pyplot.title('SSI')
#        pyplot.legend(('HS_method_1', 'GL_method_1', 'SSI_method_1', 'HS_method_2', 'GL_method_2', 'SSI_method_2'), 'best')
#        pyplot.axis([0, 2500, -2, 12])
#        pyplot.savefig(path('data/test_fitting2/SSI_GL_HS/%d_HS_GL_SSI.png') % int(id_axis_i))
        
        TT_sen_phytomer_list += TT_sen_phytomer_i_list
    
    return TT_sen_phytomer_list


def _create_absolute_TT_del_phytomer_list(id_phen_list, absolute_TT_sen_phytomer_list, second_parameters_table_dataframe):
    '''
    Create list of TT_del_phytomer.
    :Parameters:
        - `id_phen_list` : The id_phen list.
        - `absolute_TT_sen_phytomer_list` : The TT_sen_phytomer list.
        - `second_parameters_table_dataframe` : the fitted observations table, with the following headers:
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
          The table is completely filled and ordered by frequency. 
    :Types:
        - `id_phen_list` : list.
        - `absolute_TT_sen_phytomer_list` : list
        - `second_parameters_table_dataframe` : pandas.DataFrame
    '''
    TT_del_phytomer_series = pandas.Series(np.nan, index=range(len(id_phen_list)))
    df = pandas.DataFrame({'id_phen': id_phen_list, 'TT_sen_phytomer': absolute_TT_sen_phytomer_list})
    for name, group in df.groupby('id_phen'):
        second_parameters_table_dataframe_i = second_parameters_table_dataframe[second_parameters_table_dataframe['id_axis'] == name]
        a_cohort_i = second_parameters_table_dataframe_i['a_cohort'][second_parameters_table_dataframe_i.first_valid_index()]
        Nff_i = second_parameters_table_dataframe_i['Nff'][second_parameters_table_dataframe_i.first_valid_index()]
        N_cohort_i = second_parameters_table_dataframe_i['N_cohort'][second_parameters_table_dataframe_i.first_valid_index()]
        TT_del_Fhaut_phytomer_index = Nff_i - Nbr_Fhaut_persistant + N_cohort_i
        TT_del_phytomer_series[group.index[:TT_del_Fhaut_phytomer_index]] = group[:TT_del_Fhaut_phytomer_index]['TT_sen_phytomer'] + delais_phyll_sen_disp / a_cohort_i
        TT_del_phytomer_series[group.index[TT_del_Fhaut_phytomer_index:]] = TT_del_Fhaut
    return TT_del_phytomer_series.tolist()
    
    
def create_first_leaf_phen_table_dataframe(absolute_phen_table_dataframe):
    '''
    Create the first leaf phen table dataframe. 
    :Parameters:
        - `absolute_phen_table_dataframe` : the absolute phen table dataframe.
    :Types:
        - `absolute_phen_table_dataframe` : pandas.Dataframe
        
    :return: the first leaf phen table dataframe.
    :rtype: pandas.Dataframe
    '''
    assert absolute_phen_table_dataframe.count().max() == absolute_phen_table_dataframe.count().min() == absolute_phen_table_dataframe.index.size
    # Create a dataframe for first leaf (i.e. index_phytomer == 1) from absolute_phen_table_dataframe
    def first_leaf_criterion(index_i):
        # Permits to select first leaves row (i.e. row for which index_phytomer == 1).
        if absolute_phen_table_dataframe['index_phytomer'][index_i] == 1:
            return True
        return False
    first_leaf_phen_table_dataframe = absolute_phen_table_dataframe.select(first_leaf_criterion)
    return first_leaf_phen_table_dataframe


def create_phen_table_relative_dataframe(absolute_phen_table_dataframe, first_leaf_phen_table_dataframe):
    '''
    Create the relative phen table table dataframe from the absolute one. 
    :Parameters:
        - `absolute_phen_table_dataframe` : the absolute phen table dataframe.
        - `first_leaf_phen_table_dataframe` : the first leaf phen table dataframe
    :Types:
        - `absolute_phen_table_dataframe` : pandas.Dataframe
        - `first_leaf_phen_table_dataframe` : pandas.Dataframe
        
    :return: the relative phen table dataframe.
    :rtype: pandas.Dataframe
    '''
    assert absolute_phen_table_dataframe.count().max() == absolute_phen_table_dataframe.count().min() == absolute_phen_table_dataframe.index.size
    relative_phen_table_dataframe = pandas.DataFrame(index=absolute_phen_table_dataframe.index, columns=['id_phen', 'index_rel_phytomer', 'dTT_em_phytomer', 'dTT_col_phytomer', 'dTT_sen_phytomer', 'dTT_del_phytomer'])
    relative_phen_table_dataframe['id_phen'] = absolute_phen_table_dataframe['id_phen']
    tmp_series = pandas.Series(relative_phen_table_dataframe.index)
    for name, group in absolute_phen_table_dataframe.groupby('id_phen'):
        current_first_leaf_row = first_leaf_phen_table_dataframe[first_leaf_phen_table_dataframe['id_phen'] == name]
        def normalize_index_phytomer(i):
            return group['index_phytomer'][i] / float(group['index_phytomer'][group.index[-1]])
        relative_phen_table_dataframe['index_rel_phytomer'][group.index] = tmp_series[group.index].map(normalize_index_phytomer)
        def get_relative_TT_col_phytomer(i):
            return group['TT_col_phytomer'][i] - current_first_leaf_row['TT_col_phytomer'][current_first_leaf_row.first_valid_index()]
        relative_phen_table_dataframe['dTT_col_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_col_phytomer)
        def get_relative_TT_em_phytomer(i):
            return group['TT_em_phytomer'][i] - current_first_leaf_row['TT_em_phytomer'][current_first_leaf_row.first_valid_index()]
        relative_phen_table_dataframe['dTT_em_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_em_phytomer)
        def get_relative_TT_sen_phytomer(i):
            return group['TT_sen_phytomer'][i] - current_first_leaf_row['TT_sen_phytomer'][current_first_leaf_row.first_valid_index()]
        relative_phen_table_dataframe['dTT_sen_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_sen_phytomer)
        def get_relative_TT_del_phytomer(i):
            return group['TT_del_phytomer'][i] - current_first_leaf_row['TT_del_phytomer'][current_first_leaf_row.first_valid_index()]
        relative_phen_table_dataframe['dTT_del_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_del_phytomer)
        
    return relative_phen_table_dataframe


def _create_HS_GL_polynomial(HS_break_i, N_cohort_i, id_axis_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, n0_i, n1_i, n2_i, t0_i, t1_i, t2_i, a_i, c_i, d_i):
    main_stem = N_cohort_i == 1.0
    a2_i = (Nff_i - HS_break_i) / (TT_col_nff_i - TT_col_break_i)
    # define HS(TT)
    HS_1 = np.poly1d([a_cohort_i, - a_cohort_i * TT_col_0_i]) # j < HS_break_i
    HS_2 = np.poly1d([a2_i, - a2_i * TT_col_break_i + HS_break_i]) # j >= HS_break_i
    # define GL(TT) for all phases except TT < t0_i (because it depends on j)
    if main_stem:
        GL_2 = np.poly1d([(n1_i - n0_i) / (t1_i - t0_i), n0_i - t0_i * (n1_i - n0_i) / (t1_i - t0_i)])
        GL_3 = np.poly1d([(n2_i - n1_i) / (t2_i - t1_i), n1_i - t1_i * (n2_i - n1_i) / (t2_i - t1_i)])
    else: # tillers
        GL_2 = np.poly1d([n0_i])
        GL_3 = np.poly1d([(n2_i - n0_i) / (t2_i - t1_i), n0_i - t1_i * (n2_i - n0_i) / (t2_i - t1_i)])
    GL_4 = np.poly1d([a_i, - 3 * a_i * t2_i, 3 * a_i * t2_i**2 + c_i, - a_i * t2_i**3 - c_i * t2_i + d_i])
    return HS_1, HS_2, GL_2, GL_3, GL_4


def create_HS_GL_SSI_dynamic_dataframe(second_parameters_table_dataframe):
    
    HS_GL_SSI_dynamic_dataframe = pandas.DataFrame(columns=['id_axis', 'TT', 'HS', 'GL', 'SSI'])
    for i in second_parameters_table_dataframe.index:
        N_cohort_i, id_axis_i, frequency_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, dTT_MS_cohort_i, n0_i, n1_i, n2_i, t0_i, t1_i, t2_i, hs_t1_i, a_i, c_i, d_i, RMSE_gl = second_parameters_table_dataframe.ix[i].tolist()
        HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
        HS_1, HS_2, GL_2, GL_3, GL_4 = _create_HS_GL_polynomial(HS_break_i, N_cohort_i, id_axis_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, n0_i, n1_i, n2_i, t0_i, t1_i, t2_i, a_i, c_i, d_i)
        
        t0_i, t1_i, t2_i, TT_col_break_i = np.round([t0_i, t1_i, t2_i, TT_col_break_i]).astype(int)
        
        TT_1_1 = np.arange(0, t0_i)
        TT_1_2 = np.arange(t0_i, t0_i)
        TT_2_1 = np.arange(t0_i, t1_i)
        TT_2_2 = np.arange(t1_i, t1_i)
        TT_3_1 = np.arange(t1_i, t2_i)
        TT_3_2 = np.arange(t2_i, t2_i)
        TT_4_1 = np.arange(t2_i, TT_del_Fhaut)
        TT_4_2 = np.arange(TT_del_Fhaut, TT_del_Fhaut)
        
        if TT_col_break_i != 0.0: # bilinear mode
            if TT_col_break_i <= t0_i:
                TT_1_1 = np.arange(0, TT_col_break_i)
                TT_1_2 = np.arange(TT_col_break_i, t0_i)
            elif TT_col_break_i <= t1_i:
                TT_2_1 = np.arange(t0_i, TT_col_break_i)
                TT_2_2 = np.arange(TT_col_break_i, t1_i)
            elif TT_col_break_i <= t2_i:
                TT_3_1 = np.arange(t1_i, TT_col_break_i)
                TT_3_2 = np.arange(TT_col_break_i, t2_i)
            else:
                TT_4_1 = np.arange(t2_i, TT_col_break_i)
                TT_4_2 = np.arange(TT_col_break_i, TT_del_Fhaut)
            
        HS_1_TT_1_1 = np.clip(HS_1(TT_1_1), 0.0, Nff_i)
        HS_1_TT_2_1 = np.clip(HS_1(TT_2_1), 0.0, Nff_i)
        HS_1_TT_3_1 = np.clip(HS_1(TT_3_1), 0.0, Nff_i)
        HS_1_TT_4_1 = np.clip(HS_1(TT_4_1), 0.0, Nff_i)
        HS_1_TT_1_2 = np.clip(HS_1(TT_1_2), 0.0, Nff_i)
        HS_1_TT_2_2 = np.clip(HS_1(TT_2_2), 0.0, Nff_i)
        HS_1_TT_3_2 = np.clip(HS_1(TT_3_2), 0.0, Nff_i)
        HS_1_TT_4_2 = np.clip(HS_1(TT_4_2), 0.0, Nff_i)
        
        GL_1_TT_1_1 = np.clip(HS_1(TT_1_1), 0.0, 1000.0)
        GL_2_TT_2_1 = np.clip(GL_2(TT_2_1), 0.0, 1000.0)
        GL_3_TT_3_1 = np.clip(GL_3(TT_3_1), 0.0, 1000.0)
        GL_4_TT_4_1 = np.clip(GL_4(TT_4_1), 0.0, 1000.0)
        GL_1_TT_1_2 = np.clip(HS_2(TT_1_2), 0.0, 1000.0)
        GL_2_TT_2_2 = np.clip(GL_2(TT_2_2), 0.0, 1000.0)
        GL_3_TT_3_2 = np.clip(GL_3(TT_3_2), 0.0, 1000.0)
        GL_4_TT_4_2 = np.clip(GL_4(TT_4_2), 0.0, 1000.0)
        
        SSI_1_TT_1_1 = HS_1_TT_1_1 - GL_1_TT_1_1
        SSI_2_TT_2_1 = HS_1_TT_2_1 - GL_2_TT_2_1
        SSI_3_TT_3_1 = HS_1_TT_3_1 - GL_3_TT_3_1
        SSI_4_TT_4_1 = HS_1_TT_4_1 - GL_4_TT_4_1
        SSI_1_TT_1_2 = HS_1_TT_1_2 - GL_1_TT_1_2
        SSI_2_TT_2_2 = HS_1_TT_2_2 - GL_2_TT_2_2
        SSI_3_TT_3_2 = HS_1_TT_3_2 - GL_3_TT_3_2
        SSI_4_TT_4_2 = HS_1_TT_4_2 - GL_4_TT_4_2
        
        HS_GL_SSI_dynamic_dataframe_i = pandas.DataFrame(index=np.arange(TT_del_Fhaut), columns=HS_GL_SSI_dynamic_dataframe.columns)
        HS_GL_SSI_dynamic_dataframe_i['TT'] = pandas.Series(np.concatenate((TT_1_1, TT_1_2, TT_2_1, TT_2_2, TT_3_1, TT_3_2, TT_4_1, TT_4_2)))
        HS_GL_SSI_dynamic_dataframe_i['HS'] = pandas.Series(np.concatenate((HS_1_TT_1_1, HS_1_TT_1_2, HS_1_TT_2_1, HS_1_TT_2_2, HS_1_TT_3_1, HS_1_TT_3_2, HS_1_TT_4_1, HS_1_TT_4_2)))
        HS_GL_SSI_dynamic_dataframe_i['GL'] = pandas.Series(np.concatenate((GL_1_TT_1_1, GL_1_TT_1_2, GL_2_TT_2_1, GL_2_TT_2_2, GL_3_TT_3_1, GL_3_TT_3_2, GL_4_TT_4_1, GL_4_TT_4_2)))
        HS_GL_SSI_dynamic_dataframe_i['SSI'] = pandas.Series(np.concatenate((SSI_1_TT_1_1, SSI_1_TT_1_2, SSI_2_TT_2_1, SSI_2_TT_2_2, SSI_3_TT_3_1, SSI_3_TT_3_2, SSI_4_TT_4_1, SSI_4_TT_4_2)))
        HS_GL_SSI_dynamic_dataframe_i['id_axis'] = np.repeat(id_axis_i, HS_GL_SSI_dynamic_dataframe_i.index.size)
        HS_GL_SSI_dynamic_dataframe = HS_GL_SSI_dynamic_dataframe.append(HS_GL_SSI_dynamic_dataframe_i, ignore_index=True)
    
    return HS_GL_SSI_dynamic_dataframe


    