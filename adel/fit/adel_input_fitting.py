'''
This module is a facade for the functions of Adel input data fitting .

Created on Feb 7, 2012

@author: cchambon
'''  

import numpy
import pandas

import adel.fit.axis_table_fitting as axis_table_fitting
    
import adel.fit.phen_table_fitting as phen_table_fitting

# the coefficients of the secondary stem leaves number.
secondary_stem_leaves_number_coefficients = {'a_1': 0.9423, 'a_2': 0.555}
#the standard deviation used to calculate main stem emf_1 value.
emf_1_main_stem_standard_deviation = 30.0
#the Leaf number delay between Main Stem and the cohort.
leaf_number_delay_MS_cohort_dict = {'3': 1.6173, '4': 2.5181, '5': 3.4189, '6': 4.5576, '7': 5.8097}

def fit_adel_input_data_first(plant_number=100, 
                              cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.967, '6': 0.817, '7': 0.083, '8': 0.0, '9': 0.0, '10': 0.0}, 
                              main_stem_leaves_number_probability_distribution={'10': 0.145, '11': 0.818, '12': 0.036, '13': 0.0, '14': 0.0}, 
                              bolting_date=500, 
                              flowering_date=1000):
    '''
    Fit the axis table data and initialize the parameters for PhenTable fitting.
    :Parameters:
        - `plant_number` : the number of plants.
        - `cohort_probabilities` : the cohort probabilities.
        - `main_stem_leaves_number_probability_distribution` : the probability distribution of 
          the main stem leaves number.
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
    :Types:
        - `plant_number` : int.
        - `cohort_probabilities` : dict.
        - `main_stem_leaves_number_probability_distribution` : dict
        - `bolting_date` : int
        - `flowering_date` : int

    :return: The fitted axis table data and the initialized parameters for PhenTable fitting.
    :rtype: a dict of pandas.DataFrame
    '''    
    #Create and fit AxisTable 
    plant_ids = range(1,plant_number + 1)
    index_axis_list = axis_table_fitting.create_index_axis_list(plant_ids, cohort_probabilities)
    index_plt_list = axis_table_fitting.create_index_plt_list(plant_ids, index_axis_list)
    N_phyt_list = axis_table_fitting.create_N_phyt_list(index_axis_list, main_stem_leaves_number_probability_distribution, secondary_stem_leaves_number_coefficients)
    T_em_leaf1_list = axis_table_fitting.create_T_em_leaf1_list(index_axis_list, emf_1_main_stem_standard_deviation)
    # Remarque: avant de remplir la colonne TT_stop_axis il faut que la colonne TT_em_leaf1 soit totalement remplie (MB et Talles)
    T_stop_axis_list = axis_table_fitting.create_T_stop_axis_list(len(index_axis_list), int(len(index_axis_list)/2), T_em_leaf1_list, bolting_date, flowering_date)
    id_dim_list = axis_table_fitting.create_id_dim_list(index_axis_list, N_phyt_list)
    id_phen_list = axis_table_fitting.create_id_phen_list(index_axis_list, N_phyt_list)
    id_ear_list = axis_table_fitting.create_id_ear_list(index_plt_list)
    axis_array = numpy.array([index_plt_list, index_axis_list, N_phyt_list, T_stop_axis_list, id_dim_list, id_phen_list, id_ear_list, T_em_leaf1_list]).transpose()
    axis_table_dataframe = pandas.DataFrame(axis_array, columns=['id_plt', 'id_axis', 'N_phytomer', 'TT_stop_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_em_phytomer1'], dtype=float)
    
    #Initialize the parameters for PhenTable fitting.
    id_phen_without_duplicate_list = list(set(id_phen_list))
    N_cohort = [id_phen[:-2] for id_phen in id_phen_without_duplicate_list]
    axis_frequency_list = axis_table_fitting.create_axis_frequency_list(id_phen_list, id_phen_without_duplicate_list)
    Nff = [id_phen[-2:] for id_phen in id_phen_without_duplicate_list]
    a_cohort_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_col_0_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_HS_break_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_HS_NFF_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    dTT_MS_cohort_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    phen_table_parameter_array = numpy.array([N_cohort, id_phen_without_duplicate_list, axis_frequency_list, Nff, a_cohort_list, TT_col_0_list, TT_HS_break_list, TT_HS_NFF_list, dTT_MS_cohort_list]).transpose()
    unsorted_phen_table_parameter_dataframe = pandas.DataFrame(phen_table_parameter_array, columns=['N_cohort', 'id_axis', 'frequency', 'Nff', 'a_cohort', 'TT_col_0', 'TT_col_break', 'TT_col_nff', 'dTT_MS_cohort'], dtype=float)
    sorted_phen_table_parameter_dataframe = pandas.DataFrame(columns=['N_cohort', 'id_axis', 'frequency', 'Nff', 'a_cohort', 'TT_col_0', 'TT_col_break', 'TT_col_nff', 'dTT_MS_cohort'], dtype=float)
    for name, group in unsorted_phen_table_parameter_dataframe.groupby('N_cohort'):
        sorted_group = group.sort_index(by='frequency', ascending=False)
        sorted_phen_table_parameter_dataframe = sorted_phen_table_parameter_dataframe.append(sorted_group)

    return {'axis_table': axis_table_dataframe, 'parameters': sorted_phen_table_parameter_dataframe}


def fit_adel_input_data_second(first_axis_table_dataframe, user_phen_table_parameter_dataframe):
    '''
    Fit the parameters provided by the user, create the absolute and relative phen tables, complete the axis table fitting, 
    and create dim table.
    :Parameters:
        - `first_axis_table_dataframe` : the axis table from calculated with fit_adel_input_data_first.
        - `user_phen_table_parameter_dataframe` : parameters provided by the user for for PhenTable fitting.
    :Types:
        - `first_axis_table_dataframe` : pandas.DataFrame
        - `user_phen_table_parameter_dataframe` : pandas.DataFrame

    :return: The completed fitted_axis_table_dataframe, the absolute fitted phen table data, the relative fitted phen table data, 
    and the fitted parameters used for PhenTable fitting.  
    :rtype: a dict of pandas.DataFrame
    '''
    
    # Fit the parameters provided by the user
    fitted_parameter_dataframe = _fit_user_parameters(user_phen_table_parameter_dataframe, leaf_number_delay_MS_cohort_dict)
    
    # Create absolute PhenTable
    id_phen_list = phen_table_fitting.create_id_phen_list(fitted_parameter_dataframe)
    absolute_index_phytomer_list = phen_table_fitting.create_absolute_index_phytomer_list(id_phen_list)
    absolute_TT_col_phytomer_list = phen_table_fitting.create_absolute_TT_col_phytomer_list(fitted_parameter_dataframe)
    absolute_TT_em_phytomer_list = phen_table_fitting.create_absolute_TT_em_phytomer_list(absolute_TT_col_phytomer_list, fitted_parameter_dataframe)
    absolute_phen_table_array = numpy.array([id_phen_list, absolute_index_phytomer_list, absolute_TT_col_phytomer_list, absolute_TT_em_phytomer_list]).transpose()
    absolute_phen_table_dataframe = pandas.DataFrame(absolute_phen_table_array, columns=['id_phen', 'index_phytomer', 'TT_col_phytomer', 'TT_em_phytomer'], dtype=float)
    
    # Create a dataframe for first leaf (i.e. index_phytomer == 1) from absolute_phen_table_dataframe
    def first_leaf_criterion(index_i):
        # Permits to select first leaves row (i.e. row for which index_phytomer == 1).
        if absolute_phen_table_dataframe['index_phytomer'][index_i] == 1:
            return True
        return False
    
    first_leaf_phen_table_dataframe = absolute_phen_table_dataframe.select(first_leaf_criterion)
    first_leaf_phen_table_dataframe = pandas.DataFrame.from_records(first_leaf_phen_table_dataframe.to_records(index=False), index='id_phen')
    
    # Create relative PhenTable
    relative_phen_table_dataframe = phen_table_fitting.create_phen_table_relative_dataframe(absolute_phen_table_dataframe, first_leaf_phen_table_dataframe)
    #TODO: Complete AxisTable
    second_axis_table_dataframe = first_axis_table_dataframe.copy()
    
    #TODO: Create DimTable

    
    return {'axis_table': second_axis_table_dataframe, 'absolute_phen_table': absolute_phen_table_dataframe,
            'relative_phen_table': relative_phen_table_dataframe, 'parameters': fitted_parameter_dataframe}


def _fit_user_parameters(user_phen_table_parameter_dataframe, leaf_number_delay_MS_cohort_dict):
    '''
    Fit user observations. A minimal set of observations must be provided. 
    If the user does not provide complete observations, the missing observations are fitted, using the minimal set of 
    observations. In the case, extra observations (which do not belong to the minimal observations set) are ignored, and 
    will be replaced by fitted observations.  
    :Parameters:
        - `user_phen_table_parameter_dataframe` : the user observations table, with the following headers: 
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
        - `user_phen_table_parameter_dataframe` : pandas.DataFrame
        - `leaf_number_delay_MS_cohort_dict` : dict.
        
    :return: A (fitted) copy of user_phen_table_parameter_dataframe. 
    :rtype: pandas.DataFrame
    ''' 
    assert user_phen_table_parameter_dataframe['N_cohort'].count() == user_phen_table_parameter_dataframe['N_cohort'].size
    assert user_phen_table_parameter_dataframe['id_axis'].count() == user_phen_table_parameter_dataframe['id_axis'].size
    assert user_phen_table_parameter_dataframe['frequency'].count() == user_phen_table_parameter_dataframe['frequency'].size
    assert user_phen_table_parameter_dataframe['Nff'].count() == user_phen_table_parameter_dataframe['Nff'].size
    assert user_phen_table_parameter_dataframe.ix[0].count() == user_phen_table_parameter_dataframe.ix[0].size
    copy_dataframe = user_phen_table_parameter_dataframe.copy()
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

