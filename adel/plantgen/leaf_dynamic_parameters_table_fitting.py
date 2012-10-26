# -*- python -*-
#
#       Adel.PlantGen
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
from scipy.optimize import leastsq

from adel.plantgen import fit_config


def fit_user_leaf_dynamic_parameters_first(first_axis_table_id_phen_list):
    '''
    Initialize the leaf_dynamic_parameters table.
    :Parameters:
        - `first_axis_table_id_phen_list` : the number of plants.

    :Types:
        - `first_axis_table_id_phen_list` : int.

    :return: The initialized leaf_dynamic_parameters table.
    :rtype: pandas.DataFrame
    ''' 
    id_phen_without_duplicate_list = list(set(first_axis_table_id_phen_list))
    N_cohort = [float(str(int(id_phen))[:-2]) for id_phen in id_phen_without_duplicate_list]
    axis_frequency_list = _create_axis_frequency_list(first_axis_table_id_phen_list, id_phen_without_duplicate_list)
    Nff = [float(str(int(id_phen))[-2:]) for id_phen in id_phen_without_duplicate_list]
    a_cohort_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_col_0_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_HS_break_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_HS_NFF_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    dTT_MS_cohort_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    n0_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    n1_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    n2_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    t0_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    t1_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    hs_t1_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    a_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    c_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    RMSE_gl_list = [np.nan for i in range(len(id_phen_without_duplicate_list))]
    leaf_dynamic_parameters_table_array = np.array([N_cohort, id_phen_without_duplicate_list, axis_frequency_list, Nff, a_cohort_list, TT_col_0_list, TT_HS_break_list, TT_HS_NFF_list, dTT_MS_cohort_list, n0_list, n1_list, n2_list, t0_list, t1_list, hs_t1_list, a_list, c_list, RMSE_gl_list]).transpose()
    # sort leaf_dynamic_parameters table according N_cohort (ascending order) then frequency (descending order).
    unsorted_leaf_dynamic_parameters_table_dataframe = pandas.DataFrame(leaf_dynamic_parameters_table_array, columns=['N_cohort', 'id_axis', 'frequency', 'Nff', 'a_cohort', 'TT_col_0', 'TT_col_break', 'TT_col_nff', 'dTT_MS_cohort', 'n0', 'n1', 'n2', 't0', 't1', 'hs_t1', 'a', 'c', 'RMSE_gl'], dtype=float)
    sorted_leaf_dynamic_parameters_table_dataframe = pandas.DataFrame(columns=unsorted_leaf_dynamic_parameters_table_dataframe.columns, dtype=float)
    for name, group in unsorted_leaf_dynamic_parameters_table_dataframe.groupby('N_cohort'):
        sorted_group = group.sort_index(by='frequency', ascending=False)
        sorted_leaf_dynamic_parameters_table_dataframe = sorted_leaf_dynamic_parameters_table_dataframe.append(sorted_group)
    sorted_leaf_dynamic_parameters_table_dataframe.index = range(sorted_leaf_dynamic_parameters_table_dataframe.index.size)
    return sorted_leaf_dynamic_parameters_table_dataframe


def _create_axis_frequency_list(first_axis_table_id_phen_from_list, first_axis_table_id_phen_without_duplicate_list):
    '''
    Create a list of axis frequency.
    :Parameters:
        - `first_axis_table_id_phen_from_list` : the id_phen identifiers from AxisTable.
        - `first_axis_table_id_phen_without_duplicate_list` : the id_phen identifiers from AxisTable without duplicate.
    :Types:
        - `first_axis_table_id_phen_from_list` : list
        - `first_axis_table_id_phen_without_duplicate_list` : list
        
    :return: the list of axis frequency.
    :rtype: list
    '''
    axis_frequency_list = []
    for id_phen in first_axis_table_id_phen_without_duplicate_list:
        axis_frequency_list.append(first_axis_table_id_phen_from_list.count(id_phen))
    return axis_frequency_list


def fit_user_leaf_dynamic_parameters_second(user_parameter_table_dataframe, user_organ_dimensions_table_dataframe, GL_number, leaf_number_delay_MS_cohort_dict=fit_config.leaf_number_delay_MS_cohort_dict):
    '''
    Fit user observations. A minimal set of observations must be provided. 
    If the user does not provide complete observations, the missing observations are fitted, using the minimal set of 
    observations. In the case, extra observations (which do not belong to the minimal observations set) are ignored, and 
    will be replaced by fitted observations.  
    For each parameter below, we specify when it is mandatory.
    :Parameters:
        - `user_parameter_table_dataframe` : the user observations table, with the following headers: 
            * N_cohort: the cohort number. No unspecified observation.
            * id_axis: the identifier made from concatenation of N_cohort and Nff. No unspecified observation and no duplicate.
            * frequency: the occurrence frequency of id_axis. No unspecified observation.
            * Nff: the final leaves number of the current axis. No unspecified observation.
            * a_cohort: MANDATORY. The slope of phytomer emergence. This parameter can be either: 
                        - specified for each row. In this case, the a_cohort observations are not fitted.
                        - or specified for the first row of each cohort (i.e. the most frequent axis of each cohort). In this case, 
                          the set values are preserved, and the missing values are fitted using this hypothesis: 
                          a_cohort[i] = Nff[i] / (TT_col_nff[i] - TT_col_0[i]), with i the row number.
                        - or specified for the first row only. In this case, the following a_cohort observations are fitted 
                          using this hypothesis: a_cohort[i] = Nff[i] / (TT_col_nff[i] - TT_col_0[i]), with i the row number.
            * TT_col_0: MANDATORY. The Thermal Time for Haun Stage equal to 0. This parameter can be either: 
                        - specified for each row. In this case, the TT_col_0 observations are not fitted.
                        - specified for the first row of each cohort (i.e. the most frequent axis of each cohort). In this case, the missing 
                          TT_col_0 observations are fitted using this hypothesis:
                            - since same cohort axes have the same emergence date, for each cohort the value of the most frequent axis 
                              is propagated to the other axes.
                        - or specified for the first row only. In this case, the following TT_col_0 observations are fitted 
                          using this hypothesis: 
                            - same cohort axes have the same emergence date, thus the same TT_col_0 observation.
                            - TT_col_0[i] = TT_col_0[0] + (leaf_number_delay_MS_cohort_dict[N_cohort] / a_cohort[0]), 
                              with i the row number. 
            * TT_col_break: MANDATORY. The Thermal Time when the slope of phytomers emergence is changing. This parameter can be either: 
                            - specified for each row. In this case, the TT_col_break observations are not fitted.
                            - or specified for the first row only. In this case, the following TT_col_break observations are fitted 
                              using this hypothesis: TT_col_break[i] = TT_col_break[0], with i the row number.
                              This parameter is needed only for 'bilinear' calculation method.
            * TT_col_nff: MANDATORY. The Thermal Time when Haun Stage is equal to the total number of vegetative phytomers formed on an axis.
                          This parameter can be either: 
                            - specified for each row. In this case, the TT_col_nff observations are not fitted.
                            - or specified for the first row of each cohort (i.e. the most frequent axis of each cohort). In this case, 
                              the set values are preserved, and the missing values are fitted using this hypothesis: 
                              TT_col_nff[i] = dTT_MS_cohort[i] + TT_col_nff[0], with i the row number.
                            - or specified for the first row only. In this case, the following TT_col_nff observations are fitted 
                              using this hypothesis: TT_col_nff[i] = dTT_MS_cohort[i] + TT_col_nff[0], with i the row number. 
            * dTT_MS_cohort: MANDATORY. The thermal time delta between the main stem flowering and the current raw tiller flowering.
                             This parameter can be either: 
                                - specified for each row. In this case, the dTT_MS_cohort observations are not fitted.
                                - or specified for the first row of each cohort group. Since cohort group are ordered by descendant frequency, 
                                  the first row of each cohort group is also the most frequent axis of this cohort group. In this case, 
                                  for each group, the following dTT_MS_cohort observations are fitted using this hypothesis: 
                                  dTT_MS_cohort[N_cohort][j] = dTT_MS_cohort[N_cohort][0] + (id_axis[N_cohort][j] - id_axis[N_cohort][0]) / (4 * a_cohort[N_cohort][0]),
                                  with j the line number relative to N_cohort.
            * n0: MANDATORY. Number of green leaves (decimal) at t0. 
                  This parameter can be either:
                      - specified for each row. In this case, the n0 observations are not fitted.
                      - or specified for the first row only. In this case, the following n0 observations are fitted 
                        using this hypothesis: 
                            - for the main stem: n0[i] = n0[0] * Nff[i] / Nff[0], with i the row number. 
                            - for the tillers: n0[i] = min(n1[i], hs_t1[i]), with i the row number. 
            * n1: MANDATORY. Number of green leaves (decimal) at t1. 
                  This parameter can be either:
                      - specified for each row. In this case, the n1 observations are not fitted.
                      - or specified for the first row only. In this case, the following n1 observations are fitted 
                        using this hypothesis: 
                            - for the main stem: n1[i] = n1[0] * Nff[i] / Nff[0], with i the row number.
                            - for the tillers: n1[i] = n1[0], with i the row number.
            * n2: MANDATORY. Number of green leaves (decimal) at TT_col_nff. 
                  This parameter can be either:
                      - specified for each row. In this case, the n2 observations are not fitted.
                      - or specified for the first row only. In this case, the following n2 observations are fitted 
                        using this hypothesis: 
                            - for the main stem: n2[i] = n2[0] * Nff[i] / Nff[0], with i the row number.
                            - for the tillers: n2[i] = n2[0] * (1.0 - n2_MS_div_n2_cohort), with i the row number.
            * t0: Thermal time at the start of leaf senescence. 
                  The routines estimates t0 for each axis from main stem value given by the user.              
            * t1: Thermal time at the start of stem extension (bolting date). The routine estimates t1 for all axes, depending on the number of elongated internodes (add ref). 
                  - for the most frequent axis of this main stem: 
                      - in linear mode: t1[0] = TT_col_phytomer[decimal_elongated_internode_number] = TT_col_0[0] + decimal_elongated_internode_number / a_cohort[0]
                      - in bilinear mode:
                        HS_break[0] = a_cohort[0] * (TT_col_break[0] - TT_col_0[0])
                        a2[0] = (Nff[0] - HS_break[0]) / (TT_col_nff[0] - TT_col_break[0])
                        if decimal_elongated_internode_number < HS_break[0]:
                            t1[0] = TT_col_phytomer[decimal_elongated_internode_number] = TT_col_0[0] + decimal_elongated_internode_number / a_cohort[0]
                        else:
                            t1[0] = TT_col_phytomer[decimal_elongated_internode_number] = (decimal_elongated_internode_number - HS_break[0]) / a2[0] + TT_col_break[0]
                  - for all other tillers: t1[i] = t1[0] + dTT_MS_cohort[i]                        
                  with i the row number, and with decimal_elongated_internode_number defined by y = -0.0089 * x**2 + 0.4365 * x + decimal_elongated_internode_number, 
                  where x and y are respectively the internode lengths and the phytomer numbers of the most frequent axis from dimTable.
            * hs_t1: Haun Stage at t1. The routine estimates hs_t1 for each axis from main stem value given by the user.
                  hs_t1[i] = a_cohort[i] * (t1[i] - TT_col_0[i]), with i the row number.
            * a: Coefficient of the 3rd order term of the polynomial describing the dynamics of Green Leaf number after flowering. 
                 - For the main stem most frequent axis:
                   Given x[i] = GL_number.keys(), y[i] = GL_number.values(), b[0] = 0, 
                   c[0] = -((Nff[0] - decimal_elongated_internode_number) - (n2[0] - n1[0])) / (TT_col_nff[0] - t1[0])
                   and n2[0], a[0] is such as y[i] = a[0] * x[i]**3 + b[0] * x[i]**2 + c[0] * x[i] + n2[0], 
                   with i the row number, and with decimal_elongated_internode_number defined by y = -0.0089 * x**2 + 0.4365 * x + decimal_elongated_internode_number, 
                   where x and y are respectively the internode lengths and the phytomer numbers of the most frequent axis from dimTable.
                 - for the other axes: 
                   c[i] = c[0] * n2[i] / n2[0]
                   a[i] = a[0] * n2[i] / n2[0]
            * c: ??? 
                 - for the most frequent axis of the main stem: c[0] = -((Nff[0] - decimal_elongated_internode_number) - (n2[0] - n1[0])) / (TT_col_nff[0] - t1[0])
                 - for the other axes: c[i] = c[0] * n2[i] / n2[0]
            * RMSE_gl: the RMSE for the dynamic of green leaf number after estimation of parameter 'a'.
          The table is ordered by frequency.   
        - `user_organ_dimensions_table_dataframe` : the user dim table.
        - `GL_number` : the GL decimal number measured at several thermal time (including the senescence end).
        - `leaf_number_delay_MS_cohort_dict` : the Leaf number delay between Main Stem and the cohort.
    :Types:
        - `user_parameter_table_dataframe` : pandas.DataFrame
        - `user_organ_dimensions_table_dataframe` : pandas.DataFrame
        - `GL_number` : a dict of 2-tuples
        - `leaf_number_delay_MS_cohort_dict` : dict.
        
    :return: A fitted copy of user_parameter_table_dataframe. 
    :rtype: pandas.DataFrame
    ''' 
    assert user_parameter_table_dataframe['N_cohort'].count() == user_parameter_table_dataframe['N_cohort'].size
    assert user_parameter_table_dataframe['id_axis'].count() == user_parameter_table_dataframe['id_axis'].size
    assert user_parameter_table_dataframe['frequency'].count() == user_parameter_table_dataframe['frequency'].size
    assert user_parameter_table_dataframe['Nff'].count() == user_parameter_table_dataframe['Nff'].size
    MS_df = user_parameter_table_dataframe[user_parameter_table_dataframe['N_cohort'] == 1.0]
    most_frequent_MS_df = MS_df.ix[0:0]
    decimal_elongated_internode_number = _calculate_decimal_elongated_internode_number(user_organ_dimensions_table_dataframe)
    most_frequent_MS_df = _fit_most_frequent_MS_GL_dynamic(most_frequent_MS_df, decimal_elongated_internode_number, GL_number)
    other_MS_df = MS_df.ix[1:]
    other_MS_df = _fit_other_MS_HS_dynamic(most_frequent_MS_df, other_MS_df)
    other_MS_df = _fit_other_MS_GL_dynamic(most_frequent_MS_df, other_MS_df)
    tiller_axes_df = user_parameter_table_dataframe[user_parameter_table_dataframe['N_cohort'] != 1.0]
    grouped = tiller_axes_df.groupby('N_cohort')
    most_frequent_tiller_axes = []
    for N_cohort, group_indexes in grouped.groups.iteritems():
        most_frequent_tiller_axes.append(tiller_axes_df.ix[group_indexes[0:1]])
    most_frequent_tiller_axes_df = pandas.concat(most_frequent_tiller_axes)
    most_frequent_tiller_axes_df = _fit_most_frequent_tiller_axes_HS_dynamic(most_frequent_MS_df, most_frequent_tiller_axes_df)
    most_frequent_tiller_axes_df = _fit_most_frequent_tiller_axes_GL_dynamic(most_frequent_MS_df, most_frequent_tiller_axes_df)
    other_tiller_axes_df = tiller_axes_df.drop(most_frequent_tiller_axes_df.index)
    other_tiller_axes_df = _fit_other_tiller_axes_HS_dynamic(most_frequent_MS_df, most_frequent_tiller_axes_df, other_tiller_axes_df)
    other_tiller_axes_df = _fit_other_tiller_axes_GL_dynamic(most_frequent_MS_df, most_frequent_tiller_axes_df, other_tiller_axes_df)
    
    return pandas.concat([most_frequent_MS_df, other_MS_df, most_frequent_tiller_axes_df, other_tiller_axes_df]).sort()
    

def _fit_most_frequent_MS_GL_dynamic(most_frequent_MS_df, decimal_elongated_internode_number, GL_number):
    '''return fitted version of most_frequent_MS_df.'''
    # fit of t1
    most_frequent_MS_df = most_frequent_MS_df.copy()
    if most_frequent_MS_df['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_MS_df['t1'] = most_frequent_MS_df['TT_col_0'] + decimal_elongated_internode_number / most_frequent_MS_df['a_cohort']
    else: # bilinear mode
        HS_break_0 = most_frequent_MS_df['a_cohort'][0] * (most_frequent_MS_df['TT_col_break'][0] - most_frequent_MS_df['TT_col_0'][0])
        a2_0 = (most_frequent_MS_df['Nff'][0] - HS_break_0) / (most_frequent_MS_df['TT_col_nff'][0] - most_frequent_MS_df['TT_col_break'][0])
        if decimal_elongated_internode_number < HS_break_0:
            most_frequent_MS_df['t1'] = most_frequent_MS_df['TT_col_0'] + decimal_elongated_internode_number / most_frequent_MS_df['a_cohort']
        else:
            most_frequent_MS_df['t1'] = (decimal_elongated_internode_number - HS_break_0) / a2_0 + most_frequent_MS_df['TT_col_break']
    # calculation of hs_t1
    most_frequent_MS_df['hs_t1'] = most_frequent_MS_df['a_cohort'] * (most_frequent_MS_df['t1'] - most_frequent_MS_df['TT_col_0'])
    # fit of t0
    if most_frequent_MS_df['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_MS_df['t0'] = most_frequent_MS_df['TT_col_0'] + most_frequent_MS_df['n0'] / most_frequent_MS_df['a_cohort']
    else: # bilinear mode
        HS_break = most_frequent_MS_df['a_cohort'] * (most_frequent_MS_df['TT_col_break'] - most_frequent_MS_df['TT_col_0'])
        a2 = (most_frequent_MS_df['Nff'] - most_frequent_MS_df['HS_break']) / (most_frequent_MS_df['TT_col_nff'] - most_frequent_MS_df['TT_col_break'])
        n0_smaller_than_HS_break_indexes = most_frequent_MS_df[most_frequent_MS_df['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = most_frequent_MS_df[most_frequent_MS_df['n0'] >= HS_break].index
        most_frequent_MS_df['t0'][n0_smaller_than_HS_break_indexes] = most_frequent_MS_df['TT_col_0'][n0_smaller_than_HS_break_indexes] + most_frequent_MS_df['n0'][n0_smaller_than_HS_break_indexes] / most_frequent_MS_df['a_cohort'][n0_smaller_than_HS_break_indexes]
        most_frequent_MS_df['t0'][n0_greater_than_HS_break_indexes] = (most_frequent_MS_df['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + most_frequent_MS_df['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # Fit of the polynomial function describing post flowering dynamics of GL
    # c 
    most_frequent_MS_df['c'] = -((most_frequent_MS_df['Nff'] - decimal_elongated_internode_number) - (most_frequent_MS_df['n2'] - most_frequent_MS_df['n1'])) / (most_frequent_MS_df['TT_col_nff'] - most_frequent_MS_df['t1'])
    # a
    TT_col_nff_0 = most_frequent_MS_df['TT_col_nff'][0]
    n2_0 = most_frequent_MS_df['n2'][0]
    x_meas_array = np.array([TT_col_nff_0] + GL_number.keys()) - TT_col_nff_0
    y_meas_array = np.array([n2_0] + GL_number.values())
    b, c =  0.0, most_frequent_MS_df['c'][0]
    def residuals(p, y, x):
        a, = p
        err = y - peval(x, a)
        return err
    def peval(x, a):
        return np.poly1d([a, b, c, n2_0])(x)
    p0 = [-4.0e-9]
    p, cov, infodict, mesg, ier = leastsq(residuals, p0, args=(y_meas_array, x_meas_array), full_output=1)
    most_frequent_MS_df['a'][0] = p[0]
    # RMSE_gl
    chisq = (infodict['fvec']**2).sum()
    dof = len(x_meas_array) - 1 # dof is degrees of freedom
    rmse = np.sqrt(chisq / dof)
    most_frequent_MS_df['RMSE_gl'][0] = rmse
    return most_frequent_MS_df


def _fit_other_MS_HS_dynamic(most_frequent_MS_df, other_MS_df):
    '''return fitted version of other_MS_df.'''
    other_MS_df = other_MS_df.copy()  
    if other_MS_df['TT_col_0'].count() != other_MS_df['TT_col_0'].size:
        # TT_col_0
        other_MS_df['TT_col_0'] = most_frequent_MS_df['TT_col_0'][0]
        # TT_col_break
        other_MS_df['TT_col_break'] = most_frequent_MS_df['TT_col_break'][0]
        # dTT_MS_cohort
        other_MS_df['dTT_MS_cohort'] = most_frequent_MS_df['dTT_MS_cohort'][0] + (other_MS_df['id_axis'] - most_frequent_MS_df['id_axis'][0]) / (4 * most_frequent_MS_df['a_cohort'][0])
        # TT_col_nff
        other_MS_df['TT_col_nff'] = other_MS_df['dTT_MS_cohort'] + most_frequent_MS_df['TT_col_nff'][0]
        # a_cohort
        if most_frequent_MS_df['TT_col_break'][0] == 0.0: # linear mode
            other_MS_df['a_cohort'] = other_MS_df['Nff'] / (other_MS_df['TT_col_nff'] - other_MS_df['TT_col_0'])
        else: # bilinear mode
            other_MS_df['a_cohort'] = most_frequent_MS_df['a_cohort'][0] 
    return other_MS_df


def _fit_other_MS_GL_dynamic(most_frequent_MS_df, other_MS_df):
    '''return fitted version of other_MS_df.'''
    other_MS_df = other_MS_df.copy()
    # n1
    if other_MS_df['n1'].count() != other_MS_df['n1'].size:
        other_MS_df['n1'] = most_frequent_MS_df['n1'][0] * other_MS_df['Nff'] / most_frequent_MS_df['Nff'][0]
    # t1
    other_MS_df['t1'] = most_frequent_MS_df['t1'][0] + other_MS_df['dTT_MS_cohort']
    # hs_t1
    other_MS_df['hs_t1'] = other_MS_df['a_cohort'] * (other_MS_df['t1'] - other_MS_df['TT_col_0'])
    # n0 
    if other_MS_df['n0'].count() != other_MS_df['n0'].size:
        other_MS_df['n0'] = most_frequent_MS_df['n0'][0] * other_MS_df['Nff'] / most_frequent_MS_df['Nff'][0]
    # n2
    if other_MS_df['n2'].count() != other_MS_df['n2'].size:
        other_MS_df['n2'] = most_frequent_MS_df['n2'][0] * other_MS_df['Nff'] / most_frequent_MS_df['Nff'][0]
    # t0
    if most_frequent_MS_df['TT_col_break'][0] == 0.0: # linear mode
        other_MS_df['t0'] = other_MS_df['TT_col_0'] + other_MS_df['n0'] / other_MS_df['a_cohort']
    else: # bilinear mode
        HS_break = other_MS_df['a_cohort'] * (other_MS_df['TT_col_break'] - other_MS_df['TT_col_0'])
        a2 = (other_MS_df['Nff'] - other_MS_df['HS_break']) / (other_MS_df['TT_col_nff'] - other_MS_df['TT_col_break'])
        n0_smaller_than_HS_break_indexes = other_MS_df[other_MS_df['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = other_MS_df[other_MS_df['n0'] >= HS_break].index
        other_MS_df['t0'][n0_smaller_than_HS_break_indexes] = other_MS_df['TT_col_0'][n0_smaller_than_HS_break_indexes] + other_MS_df['n0'][n0_smaller_than_HS_break_indexes] / other_MS_df['a_cohort'][n0_smaller_than_HS_break_indexes]
        other_MS_df['t0'][n0_greater_than_HS_break_indexes] = (other_MS_df['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + other_MS_df['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # c
    other_MS_df['c'] = most_frequent_MS_df['c'][0] * other_MS_df['n2'] / most_frequent_MS_df['n2'][0]
    # a
    other_MS_df['a'] = most_frequent_MS_df['a'][0] * other_MS_df['n2'] / most_frequent_MS_df['n2'][0]
    # RMSE_gl
    other_MS_df['RMSE_gl'] = most_frequent_MS_df['RMSE_gl'][0]
    return other_MS_df


def _fit_most_frequent_tiller_axes_HS_dynamic(most_frequent_MS_df, most_frequent_tiller_axes_df):
    '''return fitted version of most_frequent_tiller_axes_df.'''
    most_frequent_tiller_axes_df = most_frequent_tiller_axes_df.copy()
    # TT_col_break
    most_frequent_tiller_axes_df['TT_col_break'] = most_frequent_MS_df['TT_col_break'][0]
    without_nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes_df.dropna(subset=['TT_col_0']).index
    nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes_df.index - without_nan_most_frequent_tiller_axis_indexes
    # TT_col_0
    cohorts = most_frequent_tiller_axes_df['N_cohort'].ix[nan_most_frequent_tiller_axis_indexes].astype(int).values
    leaf_number_delay_MS_cohorts = np.array([fit_config.leaf_number_delay_MS_cohort_dict[cohort] for cohort in cohorts])
    most_frequent_tiller_axes_df['TT_col_0'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_MS_df['TT_col_0'][0] + (leaf_number_delay_MS_cohorts / most_frequent_MS_df['a_cohort'][0])
    # dTT_MS_cohort: nothing to do.
    # TT_col_nff
    most_frequent_tiller_axes_df['TT_col_nff'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_tiller_axes_df['dTT_MS_cohort'].ix[nan_most_frequent_tiller_axis_indexes] + most_frequent_MS_df['TT_col_nff'][0]
    # a_cohort
    if most_frequent_MS_df['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_tiller_axes_df['a_cohort'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_tiller_axes_df['Nff'].ix[nan_most_frequent_tiller_axis_indexes] / (most_frequent_tiller_axes_df['TT_col_nff'].ix[nan_most_frequent_tiller_axis_indexes] - most_frequent_tiller_axes_df['TT_col_0'].ix[nan_most_frequent_tiller_axis_indexes])
    else: # bilinear mode
        most_frequent_tiller_axes_df['a_cohort'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_MS_df['a_cohort'][0]
    return most_frequent_tiller_axes_df


def _fit_most_frequent_tiller_axes_GL_dynamic(most_frequent_MS_df, most_frequent_tiller_axes_df):
    '''return fitted version of most_frequent_tiller_axes_df.'''
    most_frequent_tiller_axes_df = most_frequent_tiller_axes_df.copy()
    without_nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes_df.dropna(subset=['n1']).index
    nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes_df.index - without_nan_most_frequent_tiller_axis_indexes
    # n1
    most_frequent_tiller_axes_df['n1'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_MS_df['n1'][0]
    # t1
    most_frequent_tiller_axes_df['t1'] = most_frequent_MS_df['t1'][0] + most_frequent_tiller_axes_df['dTT_MS_cohort']
    # hs_t1
    most_frequent_tiller_axes_df['hs_t1'] = most_frequent_tiller_axes_df['a_cohort'] * (most_frequent_tiller_axes_df['t1'] - most_frequent_tiller_axes_df['TT_col_0'])
    # n0
    n1_hs_t1_df = most_frequent_tiller_axes_df[['n1', 'hs_t1']].ix[nan_most_frequent_tiller_axis_indexes]
    if n1_hs_t1_df.index.size != 0:
        most_frequent_tiller_axes_df['n0'].ix[nan_most_frequent_tiller_axis_indexes] = n1_hs_t1_df.apply(np.min, 1)                      
    # n2
    most_frequent_tiller_axes_df['n2'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_MS_df['n2'][0] * (1.0 - fit_config.n2_MS_div_n2_cohort)
    # t0
    if most_frequent_MS_df['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_tiller_axes_df['t0'] = most_frequent_tiller_axes_df['TT_col_0'] + most_frequent_tiller_axes_df['n0'] / most_frequent_tiller_axes_df['a_cohort']
    else: # bilinear mode
        HS_break = most_frequent_tiller_axes_df['a_cohort'] * (most_frequent_tiller_axes_df['TT_col_break'] - most_frequent_tiller_axes_df['TT_col_0'])
        a2 = (most_frequent_tiller_axes_df['Nff'] - most_frequent_tiller_axes_df['HS_break']) / (most_frequent_tiller_axes_df['TT_col_nff'] - most_frequent_tiller_axes_df['TT_col_break'])
        n0_smaller_than_HS_break_indexes = most_frequent_tiller_axes_df[most_frequent_tiller_axes_df['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = most_frequent_tiller_axes_df[most_frequent_tiller_axes_df['n0'] >= HS_break].index
        most_frequent_tiller_axes_df['t0'][n0_smaller_than_HS_break_indexes] = most_frequent_tiller_axes_df['TT_col_0'][n0_smaller_than_HS_break_indexes] + most_frequent_tiller_axes_df['n0'][n0_smaller_than_HS_break_indexes] / most_frequent_tiller_axes_df['a_cohort'][n0_smaller_than_HS_break_indexes]
        most_frequent_tiller_axes_df['t0'][n0_greater_than_HS_break_indexes] = (most_frequent_tiller_axes_df['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + most_frequent_tiller_axes_df['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # c
    most_frequent_tiller_axes_df['c'] = most_frequent_MS_df['c'][0] * most_frequent_tiller_axes_df['n2'] / most_frequent_MS_df['n2'][0]
    # a
    most_frequent_tiller_axes_df['a'] = most_frequent_MS_df['a'][0] * most_frequent_tiller_axes_df['n2'] / most_frequent_MS_df['n2'][0]
    # RMSE_gl
    most_frequent_tiller_axes_df['RMSE_gl'] = most_frequent_MS_df['RMSE_gl'][0]
    return most_frequent_tiller_axes_df
    

def _fit_other_tiller_axes_HS_dynamic(most_frequent_MS_df, most_frequent_tiller_axes_df, other_tiller_axes_df):
    '''return fitted version of other_tiller_axes_df.'''
    other_tiller_axes_df = other_tiller_axes_df.copy()
    # TT_col_break
    other_tiller_axes_df['TT_col_break'] = most_frequent_MS_df['TT_col_break'][0]    
    without_nan_other_tiller_axis_indexes = other_tiller_axes_df.dropna(subset=['TT_col_0']).index
    nan_other_tiller_axis_indexes = other_tiller_axes_df.index - without_nan_other_tiller_axis_indexes
    for name, group in other_tiller_axes_df.ix[nan_other_tiller_axis_indexes].groupby('N_cohort'):
        most_frequent_tiller_axis_idx = most_frequent_tiller_axes_df.ix[most_frequent_tiller_axes_df['N_cohort'] == name].first_valid_index()
        # TT_col_0
        other_tiller_axes_df['TT_col_0'].ix[group.index] = most_frequent_tiller_axes_df['TT_col_0'][most_frequent_tiller_axis_idx]
        # dTT_MS_cohort
        other_tiller_axes_df['dTT_MS_cohort'].ix[group.index] = \
            most_frequent_tiller_axes_df['dTT_MS_cohort'][most_frequent_tiller_axis_idx] \
            + (group['id_axis'] - most_frequent_tiller_axes_df['id_axis'][most_frequent_tiller_axis_idx]) \
            / (4 * most_frequent_MS_df['a_cohort'][0])
        if most_frequent_MS_df['TT_col_break'][0] != 0.0:
            # a_cohort, bilinear mode
            other_tiller_axes_df['a_cohort'].ix[group.index] = most_frequent_tiller_axes_df['a_cohort'][most_frequent_tiller_axis_idx]
    # TT_col_nff
    other_tiller_axes_df['TT_col_nff'].ix[nan_other_tiller_axis_indexes] = other_tiller_axes_df['dTT_MS_cohort'].ix[nan_other_tiller_axis_indexes] + most_frequent_MS_df['TT_col_nff'][0]
    # a_cohort, linear mode
    if most_frequent_MS_df['TT_col_break'][0] == 0.0: # linear mode
        other_tiller_axes_df['a_cohort'].ix[nan_other_tiller_axis_indexes] = other_tiller_axes_df['Nff'].ix[nan_other_tiller_axis_indexes] / (other_tiller_axes_df['TT_col_nff'].ix[nan_other_tiller_axis_indexes] - other_tiller_axes_df['TT_col_0'])
    return other_tiller_axes_df
    

def _fit_other_tiller_axes_GL_dynamic(most_frequent_MS_df, most_frequent_tiller_axes_df, other_tiller_axes_df):
    '''return fitted version of other_tiller_axes_df.'''
    other_tiller_axes_df = other_tiller_axes_df.copy()
    without_nan_other_tiller_axis_indexes = other_tiller_axes_df.dropna(subset=['n1']).index
    nan_other_tiller_axis_indexes = other_tiller_axes_df.index - without_nan_other_tiller_axis_indexes
    for name, group in other_tiller_axes_df.ix[nan_other_tiller_axis_indexes].groupby('N_cohort'):
        most_frequent_tiller_axis_idx = most_frequent_tiller_axes_df.ix[most_frequent_tiller_axes_df['N_cohort'] == name].first_valid_index()
        # n1
        other_tiller_axes_df['n1'].ix[group.index] = most_frequent_tiller_axes_df['n1'][most_frequent_tiller_axis_idx]
        # n2
        other_tiller_axes_df['n2'].ix[group.index] = most_frequent_tiller_axes_df['n2'][most_frequent_tiller_axis_idx]
    # t1
    other_tiller_axes_df['t1'] = most_frequent_MS_df['t1'][0] + other_tiller_axes_df['dTT_MS_cohort']
    # hs_t1
    other_tiller_axes_df['hs_t1'] = other_tiller_axes_df['a_cohort'] * (other_tiller_axes_df['t1'] - other_tiller_axes_df['TT_col_0'])
    # n0
    n1_hs_t1_df = other_tiller_axes_df[['n1', 'hs_t1']].ix[nan_other_tiller_axis_indexes]
    if n1_hs_t1_df.index.size != 0:
        other_tiller_axes_df['n0'].ix[nan_other_tiller_axis_indexes] = n1_hs_t1_df.apply(np.min, 1)                      
    # t0
    if most_frequent_MS_df['TT_col_break'][0] == 0.0: # linear mode
        other_tiller_axes_df['t0'] = other_tiller_axes_df['TT_col_0'] + other_tiller_axes_df['n0'] / other_tiller_axes_df['a_cohort']
    else: # bilinear mode
        HS_break = other_tiller_axes_df['a_cohort'] * (other_tiller_axes_df['TT_col_break'] - other_tiller_axes_df['TT_col_0'])
        a2 = (other_tiller_axes_df['Nff'] - other_tiller_axes_df['HS_break']) / (other_tiller_axes_df['TT_col_nff'] - other_tiller_axes_df['TT_col_break'])
        n0_smaller_than_HS_break_indexes = other_tiller_axes_df[other_tiller_axes_df['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = other_tiller_axes_df[other_tiller_axes_df['n0'] >= HS_break].index
        other_tiller_axes_df['t0'][n0_smaller_than_HS_break_indexes] = other_tiller_axes_df['TT_col_0'][n0_smaller_than_HS_break_indexes] + other_tiller_axes_df['n0'][n0_smaller_than_HS_break_indexes] / other_tiller_axes_df['a_cohort'][n0_smaller_than_HS_break_indexes]
        other_tiller_axes_df['t0'][n0_greater_than_HS_break_indexes] = (other_tiller_axes_df['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + other_tiller_axes_df['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # c
    other_tiller_axes_df['c'] = most_frequent_MS_df['c'][0] * other_tiller_axes_df['n2'] / most_frequent_MS_df['n2'][0]
    # a
    other_tiller_axes_df['a'] = most_frequent_MS_df['a'][0] * other_tiller_axes_df['n2'] / most_frequent_MS_df['n2'][0]
    # RMSE_gl
    other_tiller_axes_df['RMSE_gl'] = most_frequent_MS_df['RMSE_gl'][0]  
    return other_tiller_axes_df
    

def _calculate_decimal_elongated_internode_number(user_organ_dimensions_table_dataframe):
    # Calculate decimal_elongated_internode_number.
    first_axis_rows_number = int(str(int(user_organ_dimensions_table_dataframe['id_dim'][0]))[-2:])
    first_axis_L_internode_series = user_organ_dimensions_table_dataframe['L_internode'][:first_axis_rows_number]
    first_axis_L_internode_series = first_axis_L_internode_series[first_axis_L_internode_series != 0.0]
    first_axis_index_phytomer_series = user_organ_dimensions_table_dataframe['index_phytomer'][first_axis_L_internode_series.index]
    return np.polyfit(first_axis_L_internode_series, first_axis_index_phytomer_series, 2)[2]
    
    
