'''
This module is a facade for the second step of Adel input data fitting .

Created on Feb 7, 2012

@author: cchambon
'''  
import pandas

import adel.fit.axis_table_fitting as axis_table_fitting
import adel.fit.phen_table_fitting as phen_table_fitting
import adel.fit.dim_table_fitting as dim_table_fitting

#the Leaf number delay between Main Stem and the cohort.
leaf_number_delay_MS_cohort_dict = {'3': 1.6173, '4': 2.5181, '5': 3.4189, '6': 4.5576, '7': 5.8097}

def fit_adel_input_data_second(first_axis_table_dataframe, user_dim_table_dataframe, user_parameter_table_dataframe):
    '''
    Complete the axis table, the dim table and the parameters table, and create the absolute/relative phen tables.
    
    :Parameters:
        - `first_axis_table_dataframe` : the first axis table.
        - `user_dim_table_dataframe` : the user dim table.
        - `user_parameter_table_dataframe` : the user parameters table.
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
        - `first_axis_table_dataframe` : pandas.DataFrame
        - `user_dim_table_dataframe` : pandas.DataFrame
        - `user_parameter_table_dataframe` : pandas.DataFrame

    :return: The fitted axis table, the fitted absolute phen table, the fitted relative phen table,
    the fitted dim table and the fitted parameters table.
    :rtype: a tuple of pandas.DataFrame
    '''
    # Fit the parameters provided by the user
    second_parameters_table_dataframe = fit_user_parameters_second(user_parameter_table_dataframe, leaf_number_delay_MS_cohort_dict)
    # Fit absolute PhenTable
    absolute_second_phen_table_dataframe = phen_table_fitting.fit_phen_table_second(second_parameters_table_dataframe)
    # Fit relative PhenTable
    relative_second_phen_table_dataframe = phen_table_fitting.create_phen_table_relative_dataframe(absolute_second_phen_table_dataframe)
    # Fit AxisTable
    second_axis_table_dataframe = axis_table_fitting.fit_axis_table_second(first_axis_table_dataframe)
    # Fit DimTable
    dim_table_dataframe = dim_table_fitting.fit_dim_table_second(user_dim_table_dataframe, absolute_second_phen_table_dataframe)

    return second_axis_table_dataframe, absolute_second_phen_table_dataframe, relative_second_phen_table_dataframe, dim_table_dataframe, second_parameters_table_dataframe
    
    
def fit_user_parameters_second(user_parameter_table_dataframe, leaf_number_delay_MS_cohort_dict=leaf_number_delay_MS_cohort_dict):
    '''
    Fit user observations. A minimal set of observations must be provided. 
    If the user does not provide complete observations, the missing observations are fitted, using the minimal set of 
    observations. In the case, extra observations (which do not belong to the minimal observations set) are ignored, and 
    will be replaced by fitted observations.  
    :Parameters:
        - `user_parameter_table_dataframe` : the user observations table, with the following headers: 
            * N_cohort: the cohort number. No unspecified observation.
            * id_axis: the identifier made from concatenation of N_cohort and Nff. No unspecified observation and no duplicate.
            * frequency: the occurrence frequency of id_axis. No unspecified observation.
            * Nff: the final leaves number of the current axis. No unspecified observation.
            * a_cohort: The slope of phytomer emergence. This parameter can be either: 
                        - specified for each row. In this case, the a_cohort observations are not fitted.
                        - or specified for the first row only. In this case, the following a_cohort observations are fitted 
                          using this hypothesis: a_cohort[i] = Nff[i] / (TT_col_nff[i] - TT_col_0[i]), with i the row number. 
            * TT_col_0: The Thermal Time for Haun Stage equal to 0. This parameter can be either: 
                        - specified for each row. In this case, the TT_col_0 observations are not fitted.
                        - or specified for the first row only. In this case, the following TT_col_0 observations are fitted 
                          using this hypothesis: 
                            - same cohort axes have the same emergence date, thus the same TT_col_0 observation.
                            - TT_col_0[i] = TT_col_0[0] + (leaf_number_delay_MS_cohort_dict[N_cohort] / a_cohort[0]), 
                              with i the row number. 
            * TT_col_break: The Thermal Time when the slope of phytomers emergence is changing. This parameter can be either: 
                            - specified for each row. In this case, the TT_col_break observations are not fitted.
                            - or specified for the first row only. In this case, the following TT_col_break observations are fitted 
                              using this hypothesis: TT_col_break[i] = TT_col_break[0], with i the row number.
                              This parameter is needed only for 'bilinear' calculation method.
            * TT_col_nff: The Thermal Time when Haun Stage is equal to the total number of vegetative phytomers formed on an axis.
                          This parameter can be either: 
                            - specified for each row. In this case, the TT_col_nff observations are not fitted.
                            - or specified for the first row only. In this case, the following TT_col_nff observations are fitted 
                              using this hypothesis: TT_col_nff[i] = dTT_MS_cohort[i] + TT_col_nff[i], with i the row number. 
            * dTT_MS_cohort: The thermal time delta between the main stem flowering and the current raw tiller flowering.
                             This parameter can be either: 
                                - specified for each row. In this case, the dTT_MS_cohort observations are not fitted.
                                - or specified for the first row of each cohort group. Since cohort group are ordered by descendant frequency, 
                                  the first row of each cohort group is also the most frequent axis of this cohort group. In this case, 
                                  for each group, the following dTT_MS_cohort observations are fitted using this hypothesis: 
                                  dTT_MS_cohort[N_cohort][j] = dTT_MS_cohort[N_cohort][0] + (id_axis[N_cohort][j] - id_axis[N_cohort][0]) / (4 * a_cohort[N_cohort][0]),
                                  with j the line number relative to N_cohort.
          The table is ordered by frequency.   
        - `leaf_number_delay_MS_cohort_dict` : the Leaf number delay between Main Stem and the cohort.
    :Types:
        - `user_parameter_table_dataframe` : pandas.DataFrame
        - `leaf_number_delay_MS_cohort_dict` : dict.
        
    :return: A fitted copy of user_parameter_table_dataframe. 
    :rtype: pandas.DataFrame
    ''' 
    assert user_parameter_table_dataframe['N_cohort'].count() == user_parameter_table_dataframe['N_cohort'].size
    assert user_parameter_table_dataframe['id_axis'].count() == user_parameter_table_dataframe['id_axis'].size
    assert user_parameter_table_dataframe['frequency'].count() == user_parameter_table_dataframe['frequency'].size
    assert user_parameter_table_dataframe['Nff'].count() == user_parameter_table_dataframe['Nff'].size
    assert user_parameter_table_dataframe.ix[0].count() == user_parameter_table_dataframe.ix[0].size
    copy_dataframe = user_parameter_table_dataframe.copy()
    tmp_series = pandas.Series(copy_dataframe.index)
    if copy_dataframe['TT_col_0'].count() != copy_dataframe['TT_col_0'].size:
        # Fit copy_dataframe['TT_col_0'].
        for name, group in copy_dataframe.groupby('N_cohort'):
            if name != 1.0:
                copy_dataframe['TT_col_0'][group.index[0]] = copy_dataframe['TT_col_0'][0] + (leaf_number_delay_MS_cohort_dict[str(int(group['N_cohort'][0]))] / copy_dataframe['a_cohort'][0])
        copy_dataframe['TT_col_0'] = copy_dataframe['TT_col_0'].fillna()
        
    if copy_dataframe['TT_col_break'].count() != copy_dataframe['TT_col_break'].size:
        # Fit copy_dataframe['TT_col_break'].
        copy_dataframe['TT_col_break'] = copy_dataframe['TT_col_break'].fillna()
        
    if copy_dataframe['dTT_MS_cohort'].count() != copy_dataframe['dTT_MS_cohort'].size:
        # Fit copy_dataframe['dTT_MS_cohort'].
        for name, group in copy_dataframe.groupby('N_cohort'):
            def fit_dTT_MS_cohort(i):
                if i == group['dTT_MS_cohort'].index[0]:
                    result = group['dTT_MS_cohort'][0]
                else:
                    result =  group['dTT_MS_cohort'][0] + (group['id_axis'][i] - group['id_axis'][0]) / (4 * copy_dataframe['a_cohort'][0])
                return result
            copy_dataframe['dTT_MS_cohort'][group.index] = tmp_series[group.index].map(fit_dTT_MS_cohort)
        
    if copy_dataframe['TT_col_nff'].count() != copy_dataframe['TT_col_nff'].size:
        # Fit copy_dataframe['TT_col_nff'].
        def fit_TT_col_nff(i):
            if i == 0:
                result = copy_dataframe['TT_col_nff'][0]
            else:
                result = copy_dataframe['dTT_MS_cohort'][i] + copy_dataframe['TT_col_nff'][0]
            return result
        copy_dataframe['TT_col_nff'][tmp_series.index] = tmp_series.map(fit_TT_col_nff)
        
    if copy_dataframe['a_cohort'].count() != copy_dataframe['a_cohort'].size:
        # Fit copy_dataframe['a_cohort'].
        if copy_dataframe['TT_col_break'][0] == 0.0: # linear mode
            def fit_a_cohort(i):
                if i == 0:
                    result = copy_dataframe['a_cohort'][0]
                else:
                    result = copy_dataframe['Nff'][i] / (copy_dataframe['TT_col_nff'][i] - copy_dataframe['TT_col_0'][i])
                return result
            copy_dataframe['a_cohort'][tmp_series.index] = tmp_series.map(fit_a_cohort)
        else: # bilinear mode
            copy_dataframe['a_cohort'] = copy_dataframe['a_cohort'].fillna()
        
    return copy_dataframe

