'''
This module is a facade for the second step of Adel input data fitting .

Created on Feb 7, 2012

@author: cchambon
'''  
import pandas
import numpy as np
from scipy.optimize import leastsq

import adel.fit.axis_table_fitting as axis_table_fitting
import adel.fit.phen_table_fitting as phen_table_fitting
import adel.fit.dim_table_fitting as dim_table_fitting

#the Leaf number delay between Main Stem and the cohort.
leaf_number_delay_MS_cohort_dict = {'3': 1.6173, '4': 2.5181, '5': 3.4189, '6': 4.5576, '7': 5.8097}
n2_MS_div_n2_cohort = 0.15


def fit_adel_input_data_second(first_axis_table_dataframe, user_dim_table_dataframe, user_parameter_table_dataframe, 
                               GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                               bolting_date=500, flowering_date=1440, 
                               delais_TT_stop_del_axis=600,
                               final_axes_number=250):
    '''
    Complete the axis table, the dim table and the parameters table, and create the absolute/relative phen tables.
    Construct:
        - the parameter table which contains the information about tillers behaviour. The parameters table is not an input of ADEL.
          It is used only to build the other tables: ParametersTable
        - ADEL input data tables: second_axis_table_dataframe, relative_second_phen_table_dataframe, relative_dim_table_dataframe
        - tables in order the user could check the fitting results: first_leaf_phen_table_dataframe, 
          HS_GL_SSI_dynamic_dataframe
    
    :Parameters:
        - `first_axis_table_dataframe` : the first axis table.
        - `user_dim_table_dataframe` : the user dim table.
        - `user_parameter_table_dataframe` : the user parameters table.
        - `GL_number` : the GL decimal number measured at several thermal time (including the senescence end).
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
        - `delais_TT_stop_del_axis` : Thermal time during which a tiller remains present on the plant after the tiller has stopped 
                                      growing (for tillers that senesce).
        - `final_axes_number` : the final number of axes per square meter.
    :Types:
        - `first_axis_table_dataframe` : pandas.DataFrame
        - `user_dim_table_dataframe` : pandas.DataFrame
        - `user_parameter_table_dataframe` : pandas.DataFrame
        - `GL_number` : a dict of 2-tuples
        - `bolting_date` : int
        - `flowering_date` : int
        - `delais_TT_stop_del_axis` : int
        - `final_axes_number` : int

    :return: The fitted axis table, the fitted absolute phen table, the fitted relative phen table,
    the fitted dim table and the fitted parameters table.
    :rtype: a tuple of pandas.DataFrame
    '''
    # Fit the parameters provided by the user
    second_parameters_table_dataframe = fit_user_parameters_second(user_parameter_table_dataframe, user_dim_table_dataframe, GL_number, leaf_number_delay_MS_cohort_dict)
    # Fit absolute PhenTable
    absolute_second_phen_table_dataframe = phen_table_fitting.fit_phen_table_second(second_parameters_table_dataframe)
    # Extract the first leaves data from absolute_second_phen_table_dataframe
    first_leaf_phen_table_dataframe = phen_table_fitting.create_first_leaf_phen_table_dataframe(absolute_second_phen_table_dataframe)
    # Fit relative PhenTable
    relative_second_phen_table_dataframe = phen_table_fitting.create_phen_table_relative_dataframe(absolute_second_phen_table_dataframe, first_leaf_phen_table_dataframe)
    # Fit AxisTable
    second_axis_table_dataframe = axis_table_fitting.fit_axis_table_second(first_axis_table_dataframe, first_leaf_phen_table_dataframe, bolting_date, flowering_date, delais_TT_stop_del_axis, final_axes_number)
    # Fit DimTable
    absolute_dim_table_dataframe = dim_table_fitting.fit_dim_table_second(user_dim_table_dataframe, absolute_second_phen_table_dataframe)
    # Fit relative dimTable
    relative_dim_table_dataframe = dim_table_fitting.create_dim_table_relative_dataframe(absolute_dim_table_dataframe)
    # Create a table with the following columns: id_axis,TT,HS,GL,SSI
    HS_GL_SSI_dynamic_dataframe = phen_table_fitting.create_HS_GL_SSI_dynamic_dataframe(second_parameters_table_dataframe)
    
    return second_axis_table_dataframe, absolute_second_phen_table_dataframe, relative_second_phen_table_dataframe, \
           absolute_dim_table_dataframe, second_parameters_table_dataframe, first_leaf_phen_table_dataframe, \
           HS_GL_SSI_dynamic_dataframe, relative_dim_table_dataframe
    
    
def fit_user_parameters_second(user_parameter_table_dataframe, user_dim_table_dataframe, GL_number, leaf_number_delay_MS_cohort_dict=leaf_number_delay_MS_cohort_dict):
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
            * n2: MANDATORY. Number of green leaves (decimal) at t2. 
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
            * t2: Thermal time at flag leaf ligulation. The routines estimates t2 for each axis from main stem value given by the user.  
                  t2[i] = TT_col_nff[i], with i the row number.
            * hs_t1: Haun Stage at t1. The routine estimates hs_t1 for each axis from main stem value given by the user.
                  hs_t1[i] = a_cohort[i] * (t1[i] - TT_col_0[i]), with i the row number.
            * a: Coefficient of the 
                 - For the main stem most frequent axis:
                   Given x[i] = GL_number.keys(), y[i] = GL_number.values(), b[0] = 0, 
                   c[0] = -((Nff[0] - decimal_elongated_internode_number) - (n2[0] - n1[0])) / (t2[0] - t1[0])
                   and d[0] = n2[0], a[0] is such as y[i] = a[0] * x[i]**3 + b[0] * x[i]**2 + c[0] * x[i] + d[0], 
                   with i the row number, and with decimal_elongated_internode_number defined by y = -0.0089 * x**2 + 0.4365 * x + decimal_elongated_internode_number, 
                   where x and y are respectively the internode lengths and the phytomer numbers of the most frequent axis from dimTable.
                 - for the other axes: 
                   c[i] = c[0] * d[i] / d[0]
                   d[i] = n2[i]
                   a[i] = a[0] * d[i] / d[0]
            * c: ??? 
                 - for the most frequent axis of the main stem: c[0] = -((Nff[0] - decimal_elongated_internode_number) - (n2[0] - n1[0])) / (t2[0] - t1[0])
                 - for the other axes: c[i] = c[0] * d[i] / d[0]
            * d: ??? 
                This parameter don't have to be specified.
                d = n2
            * RMSE_gl: the RMSE for the dynamic of green leaf number after estimation of parameter 'a'.
          The table is ordered by frequency.   
        - `user_dim_table_dataframe` : the user dim table.
        - `GL_number` : the GL decimal number measured at several thermal time (including the senescence end).
        - `leaf_number_delay_MS_cohort_dict` : the Leaf number delay between Main Stem and the cohort.
    :Types:
        - `user_parameter_table_dataframe` : pandas.DataFrame
        - `user_dim_table_dataframe` : pandas.DataFrame
        - `GL_number` : a dict of 2-tuples
        - `leaf_number_delay_MS_cohort_dict` : dict.
        
    :return: A fitted copy of user_parameter_table_dataframe. 
    :rtype: pandas.DataFrame
    ''' 
    assert user_parameter_table_dataframe['N_cohort'].count() == user_parameter_table_dataframe['N_cohort'].size
    assert user_parameter_table_dataframe['id_axis'].count() == user_parameter_table_dataframe['id_axis'].size
    assert user_parameter_table_dataframe['frequency'].count() == user_parameter_table_dataframe['frequency'].size
    assert user_parameter_table_dataframe['Nff'].count() == user_parameter_table_dataframe['Nff'].size
    copy_dataframe = user_parameter_table_dataframe.copy()
    _fit_TT_col_0_series(copy_dataframe)
    _fit_TT_col_break_series(copy_dataframe)
    _fit_dTT_MS_cohort_series(copy_dataframe)
    _fit_TT_col_nff_series(copy_dataframe)
    _fit_a_cohort_series(copy_dataframe)
    _fit_n1_series(copy_dataframe)
    decimal_elongated_internode_number = _calculate_decimal_elongated_internode_number(user_dim_table_dataframe)
    _fit_t1_series(copy_dataframe, decimal_elongated_internode_number)
    _fit_hs_t1_series(copy_dataframe)
    _fit_n0_series(copy_dataframe)
    _fit_n2_series(copy_dataframe)
    _fit_t0_series(copy_dataframe)
    _fit_t2_series(copy_dataframe)
    _fit_d_series(copy_dataframe)
    _fit_c_series(copy_dataframe, decimal_elongated_internode_number)
    _fit_a_series(copy_dataframe, GL_number)
        
    return copy_dataframe


def _fit_TT_col_0_series(copy_dataframe):
    if copy_dataframe['TT_col_0'].count() != copy_dataframe['TT_col_0'].size:
        # Fit copy_dataframe['TT_col_0'].
        for name, group in copy_dataframe.groupby('N_cohort'):
            if name != 1.0 and np.isnan(copy_dataframe['TT_col_0'][group.index[0]]):
                copy_dataframe['TT_col_0'][group.index[0]] = copy_dataframe['TT_col_0'][0] + (leaf_number_delay_MS_cohort_dict[str(int(group['N_cohort'][group.index[0]]))] / copy_dataframe['a_cohort'][0])
        copy_dataframe['TT_col_0'] = copy_dataframe['TT_col_0'].fillna()
        
        
def _fit_TT_col_break_series(copy_dataframe):
    if copy_dataframe['TT_col_break'].count() != copy_dataframe['TT_col_break'].size:
        # Fit copy_dataframe['TT_col_break'].
        copy_dataframe['TT_col_break'] = copy_dataframe['TT_col_break'].fillna()    
    
    
def _fit_dTT_MS_cohort_series(copy_dataframe):
    if copy_dataframe['dTT_MS_cohort'].count() != copy_dataframe['dTT_MS_cohort'].size:
        tmp_series = pandas.Series(copy_dataframe.index)
        # Fit copy_dataframe['dTT_MS_cohort'].
        def fit_dTT_MS_cohort(i):
            if i == group.index[0]:
                result = group['dTT_MS_cohort'][group.index[0]]
            else:
                result =  group['dTT_MS_cohort'][group.index[0]] + (group['id_axis'][i] - group['id_axis'][group.index[0]]) / (4 * copy_dataframe['a_cohort'][0])
            return result
        for name, group in copy_dataframe.groupby('N_cohort'):
            copy_dataframe['dTT_MS_cohort'][group.index] = tmp_series[group.index].map(fit_dTT_MS_cohort)


def _fit_TT_col_nff_series(copy_dataframe):
    if copy_dataframe['TT_col_nff'].count() != copy_dataframe['TT_col_nff'].size:
        tmp_series = pandas.Series(copy_dataframe.index)
        # Fit copy_dataframe['TT_col_nff'].
        def fit_TT_col_nff(i):
            if i == group.index[0] and not np.isnan(group['TT_col_nff'][group.index[0]]):
                result = copy_dataframe['TT_col_nff'][group.index[0]]
            else:
                result = copy_dataframe['dTT_MS_cohort'][i] + copy_dataframe['TT_col_nff'][0]
            return result
        for name, group in copy_dataframe.groupby('N_cohort'):
            copy_dataframe['TT_col_nff'][group.index] = tmp_series[group.index].map(fit_TT_col_nff)


def _fit_a_cohort_series(copy_dataframe):
    if copy_dataframe['a_cohort'].count() != copy_dataframe['a_cohort'].size:
        tmp_series = pandas.Series(copy_dataframe.index)
        # Fit copy_dataframe['a_cohort'].
        if copy_dataframe['TT_col_break'][0] == 0.0: # linear mode
            def fit_a_cohort(i):
                if i == group.index[0] and not np.isnan(group['a_cohort'][group.index[0]]):
                    result = group['a_cohort'][group.index[0]]
                else:
                    result = group['Nff'][i] / (group['TT_col_nff'][i] - group['TT_col_0'][i])
                return result
            for name, group in copy_dataframe.groupby('N_cohort'):
                copy_dataframe['a_cohort'][group.index] = tmp_series[group.index].map(fit_a_cohort)
        else: # bilinear mode
            copy_dataframe['a_cohort'] = copy_dataframe['a_cohort'].fillna()


def _fit_n1_series(copy_dataframe):
    if copy_dataframe['n1'].count() != copy_dataframe['n1'].size:
        # Fit copy_dataframe['n1'].
        for name, group in copy_dataframe.groupby('N_cohort'):
            if name == 1.0: # main stem
                copy_dataframe['n1'][group.index[1:]] = copy_dataframe['n1'][0] * copy_dataframe['Nff'][group.index[1:]] / copy_dataframe['Nff'][0]
            else: # tillers
                #TODO: check this with Bruno
                if np.isnan(group['n1'][group.index[0]]):
                    copy_dataframe['n1'][group.index] = copy_dataframe['n1'][0]
                else:
                    copy_dataframe['n1'][group.index[1:]] = group['n1'][group.index[0]]

def _calculate_decimal_elongated_internode_number(user_dim_table_dataframe):
    # Calculate decimal_elongated_internode_number.
    first_axis_rows_number = int(str(int(user_dim_table_dataframe['id_dim'][0]))[-2:])
    first_axis_L_internode_series = user_dim_table_dataframe['L_internode'][:first_axis_rows_number]
    first_axis_L_internode_series = first_axis_L_internode_series[first_axis_L_internode_series != 0.0]
    first_axis_index_phytomer_series = user_dim_table_dataframe['index_phytomer'][first_axis_L_internode_series.index]
    return np.polyfit(first_axis_L_internode_series, first_axis_index_phytomer_series, 2)[2]


def _fit_t1_series(copy_dataframe, decimal_elongated_internode_number):
    # Fit copy_dataframe['t1'].
    # For the most frequent axis of this main stem
    if copy_dataframe['TT_col_break'][0] == 0.0: # linear mode
        copy_dataframe['t1'][0] = copy_dataframe['TT_col_0'][0] + decimal_elongated_internode_number / copy_dataframe['a_cohort'][0]
    else: # bilinear mode
        HS_break_0 = copy_dataframe['a_cohort'][0] * (copy_dataframe['TT_col_break'][0] - copy_dataframe['TT_col_0'][0])
        a2_0 = (copy_dataframe['Nff'][0] - HS_break_0) / (copy_dataframe['TT_col_nff'][0] - copy_dataframe['TT_col_break'][0])
        if decimal_elongated_internode_number < HS_break_0:
            copy_dataframe['t1'][0] = copy_dataframe['TT_col_0'][0] + decimal_elongated_internode_number / copy_dataframe['a_cohort'][0]
        else:
            copy_dataframe['t1'][0] = (decimal_elongated_internode_number - HS_break_0) / a2_0 + copy_dataframe['TT_col_break'][0]
    # for the other axes
    copy_dataframe['t1'][1:] = copy_dataframe['t1'][0] + copy_dataframe['dTT_MS_cohort'][1:]


def _fit_hs_t1_series(copy_dataframe):
    # Fit copy_dataframe['hs_t1'].
    copy_dataframe['hs_t1'] = copy_dataframe['a_cohort'] * (copy_dataframe['t1'] - copy_dataframe['TT_col_0'])


def _fit_n0_series(copy_dataframe):
    if copy_dataframe['n0'].count() != copy_dataframe['n0'].size:
        # Fit copy_dataframe['n0'].
        for name, group in copy_dataframe.groupby('N_cohort'):
            if name == 1.0: # main stem
                copy_dataframe['n0'][group.index[1:]] = copy_dataframe['n0'][0] * copy_dataframe['Nff'][group.index[1:]] / copy_dataframe['Nff'][0]
            else: # tillers
                #TODO: check this with Bruno
                n1_hs_t1_dataframe = group[['n1', 'hs_t1']]
                copy_dataframe['n0'][group.index] = n1_hs_t1_dataframe.apply(np.min, 1)


def _fit_n2_series(copy_dataframe):
    if copy_dataframe['n2'].count() != copy_dataframe['n2'].size:
        # Fit copy_dataframe['n2'].
        for name, group in copy_dataframe.groupby('N_cohort'):
            if name == 1.0: # main stem
                copy_dataframe['n2'][group.index[1:]] = copy_dataframe['n2'][0] * copy_dataframe['Nff'][group.index[1:]] / copy_dataframe['Nff'][0]
            else: # tillers
                #TODO: check this with Bruno
                if np.isnan(group['n2'][group.index[0]]):
                    copy_dataframe['n2'][group.index] = copy_dataframe['n2'][0] * (1.0 - n2_MS_div_n2_cohort)
                else:
                    copy_dataframe['n2'][group.index[1:]] = group['n2'][group.index[0]] * (1.0 - n2_MS_div_n2_cohort)
                    

def _fit_t0_series(copy_dataframe):
    '''
    Fits t0 according to the rate of HS progress: linear or bilinear.
    - in linear mode: t0[i] = TT_col_phytomer[n0[i]] = TT_col_0[i] + n0[i] / a_cohort[i]
    - in bilinear mode: 
        HS_break[i] = a_cohort[i] * (TT_col_break[i] - TT_col_0[i])
        a2[i] = (Nff[i] - HS_break[i]) / (TT_col_nff[i] - TT_col_break[i])
        if n0[i] < HS_break[i]:
            t0[i] = TT_col_phytomer[n0[i]] = TT_col_0[i] + n0[i] / a_cohort[i]
        else:
            t0[i] = TT_col_phytomer[n0[i]] = (n0[i] - HS_break[i]) / a2[i] + TT_col_break[i]
        with i the row number. 
    '''
    # Fit copy_dataframe['t0'].
    if copy_dataframe['TT_col_break'][0] == 0.0: # linear mode
        copy_dataframe['t0'] = copy_dataframe['TT_col_0'] + copy_dataframe['n0'] / copy_dataframe['a_cohort']
    else: # bilinear mode
        HS_break = copy_dataframe['a_cohort'] * (copy_dataframe['TT_col_break'] - copy_dataframe['TT_col_0'])
        a2 = (copy_dataframe['Nff'] - copy_dataframe['HS_break']) / (copy_dataframe['TT_col_nff'] - copy_dataframe['TT_col_break'])
        n0_smaller_than_HS_break_indexes = copy_dataframe[copy_dataframe['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = copy_dataframe[copy_dataframe['n0'] >= HS_break].index
        copy_dataframe['t0'][n0_smaller_than_HS_break_indexes] = copy_dataframe['TT_col_0'][n0_smaller_than_HS_break_indexes] + copy_dataframe['n0'][n0_smaller_than_HS_break_indexes] / copy_dataframe['a_cohort'][n0_smaller_than_HS_break_indexes]
        copy_dataframe['t0'][n0_greater_than_HS_break_indexes] = (copy_dataframe['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + copy_dataframe['TT_col_break'][n0_greater_than_HS_break_indexes]    


def _fit_t2_series(copy_dataframe):
    # Fit copy_dataframe['t2'].
    copy_dataframe['t2'] = copy_dataframe['TT_col_nff']


def _fit_d_series(copy_dataframe):
    copy_dataframe['d'] = copy_dataframe['n2']


def _fit_c_series(copy_dataframe, decimal_elongated_internode_number):
    copy_dataframe['c'][0] = -((copy_dataframe['Nff'][0] - decimal_elongated_internode_number) - (copy_dataframe['n2'][0] - copy_dataframe['n1'][0])) / (copy_dataframe['t2'][0] - copy_dataframe['t1'][0])
    copy_dataframe['c'][1:] = copy_dataframe['c'][0] * copy_dataframe['d'][1:] / copy_dataframe['d'][0]
    

def _fit_a_series(copy_dataframe, GL_number):
    t2_0 = copy_dataframe['t2'][0]
    n2_0 = copy_dataframe['n2'][0]
    x_meas_array = np.array([t2_0] + GL_number.keys()) - t2_0
    y_meas_array = np.array([n2_0] + GL_number.values())

    b, c, d =  0.0, -0.00458844157413, 5.8
    
    def residuals(p, y, x):
        a, = p
        err = y - peval(x, a)
        return err
    
    def peval(x, a):
        return np.poly1d([a, b, c, d])(x)
    
    p0 = [-4.0e-9]
    
    p, cov, infodict, mesg, ier = leastsq(residuals, p0, args=(y_meas_array, x_meas_array), full_output=1)
    chisq = (infodict['fvec']**2).sum()
    # dof is degrees of freedom
    dof = len(x_meas_array) - 1
    rmse = np.sqrt(chisq / dof)

    copy_dataframe['a'][0] = p[0]

    # for the other axes
    copy_dataframe['a'][1:] = copy_dataframe['a'][0] * copy_dataframe['d'][1:] / copy_dataframe['d'][0]
    copy_dataframe['RMSE_gl'] = rmse
    
    
