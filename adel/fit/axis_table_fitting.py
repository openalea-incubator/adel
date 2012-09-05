'''
This module provides functions to calculate AxisTable.

Created on 28 nov. 2011

@author: cchambon
'''

import random
import math
import numpy as np
import pandas

# the coefficients of the secondary stem leaves number.
secondary_stem_leaves_number_coefficients = {'a_1': 0.9423, 'a_2': 0.555}
#the standard deviation used to calculate main stem emf_1 value.
emf_1_main_stem_standard_deviation = 30.0


def generate_axes(plant_number, cohort_probabilities, main_stem_leaves_number_probability_distribution):
    '''
    Fit the axis table: first step.
    Compute the following columns: 
        - id_axis: the axis ids for each plant.   
        - id_plt: the plant ids.
        - N_phytomer: the final leave numbers for each axis.
        - id_dim: the keys referring to a group of lines in the table of the 
          organ dimensions.
        - id_phen: the keys referring to a group of lines in the phenological 
          table. 
        - id_ear: the keys referring to a group of lines in the table of ear appearance dates.    
    :Parameters:
        - `plant_number` : the number of plants.
        - `cohort_probabilities` : the cohort probabilities.
        - `main_stem_leaves_number_probability_distribution` : the probability distribution of 
          the main stem leaves number.
    :Types:
        - `plant_number` : int.
        - `cohort_probabilities` : dict.
        - `main_stem_leaves_number_probability_distribution` : dict
        
    :return: The first axis table.
    :rtype: pandas.DataFrame
    '''
    plant_ids = range(1,plant_number + 1)
    index_axis_list = _create_index_axis_list(plant_ids, cohort_probabilities)
    index_plt_list = _create_index_plt_list(plant_ids, index_axis_list)
    N_phyt_list = _create_N_phyt_list(index_axis_list, main_stem_leaves_number_probability_distribution, secondary_stem_leaves_number_coefficients)
    TT_stop_axis_list = [np.nan for i in range(len(index_axis_list))]
    TT_del_axis_list = [np.nan for i in range(len(index_axis_list))]
    id_dim_list = _create_id_dim_list(index_axis_list, N_phyt_list)
    id_phen_list = _create_id_phen_list(index_axis_list, N_phyt_list)
    id_ear_list = _create_id_ear_list(index_plt_list)
    TT_em_phytomer1_list = [np.nan for i in range(len(index_axis_list))]
    TT_col_phytomer1_list = [np.nan for i in range(len(index_axis_list))]
    TT_sen_phytomer1_list = [np.nan for i in range(len(index_axis_list))]
    TT_del_phytomer1_list = [np.nan for i in range(len(index_axis_list))]
    axis_table_array = np.array([index_plt_list, index_axis_list, N_phyt_list, TT_stop_axis_list, TT_del_axis_list, id_dim_list, id_phen_list, id_ear_list, TT_em_phytomer1_list, TT_col_phytomer1_list, TT_sen_phytomer1_list, TT_del_phytomer1_list]).transpose()
    return pandas.DataFrame(axis_table_array, columns=['id_plt', 'id_axis', 'N_phytomer', 'TT_stop_axis', 'TT_del_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_em_phytomer1', 'TT_col_phytomer1', 'TT_sen_phytomer1', 'TT_del_phytomer1'], dtype=float)


def fit_axis_table_second(first_axis_table_dataframe, first_leaf_phen_table_dataframe, bolting_date, flowering_date, delais_TT_stop_del_axis, final_axes_number):
    '''
    Fit the axis table: second step.
    :Parameters:
        - `first_axis_table_dataframe` : The first axis table.
        - `first_leaf_phen_table_dataframe` : the first leaf phen table.
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
        - `delais_TT_stop_del_axis` : ???
        - `final_axes_number` : the final number of axes
    :Types:
        - `first_axis_table_dataframe` : pandas.DataFrame
        - `first_leaf_phen_table_dataframe` : pandas.Dataframe
        - `bolting_date` : int
        - `flowering_date` : int
        - `delais_TT_stop_del_axis` : int
        - `final_axes_number` : int
        
    :return: The second axis table.
    :rtype: pandas.DataFrame
    '''
    second_axis_table_dataframe = first_axis_table_dataframe.copy()
    (second_axis_table_dataframe['TT_em_phytomer1'], 
     second_axis_table_dataframe['TT_col_phytomer1'], 
     second_axis_table_dataframe['TT_sen_phytomer1'],
     second_axis_table_dataframe['TT_del_phytomer1']) = _create_all_TT_phytomer1_list(first_axis_table_dataframe, emf_1_main_stem_standard_deviation, first_leaf_phen_table_dataframe)
    second_axis_table_dataframe['TT_stop_axis'] = _create_TT_stop_axis_list(first_axis_table_dataframe.index.size, final_axes_number, second_axis_table_dataframe['TT_em_phytomer1'].tolist(), bolting_date, flowering_date)
    second_axis_table_dataframe['TT_del_axis'] = _create_TT_del_axis_list(second_axis_table_dataframe['TT_stop_axis'], delais_TT_stop_del_axis)
    
    return second_axis_table_dataframe
    

def _create_index_plt_list(plant_ids, first_axis_table_index_axis_list):
    '''
    Create plant indexes column.
    :Parameters:
        - `plant_ids` : the plant indexes.
        - `first_axis_table_index_axis_list` : the axes column.
    :Types:
        - `plant_ids` : list
        - `first_axis_table_index_axis_list` : list
        
    :return: The plant indexes column.
    :rtype: list
    '''
    index_plt_list = []
    current_plant_index = 0
    for plant_id in plant_ids:
        start_index = current_plant_index + 1
        if 1 in first_axis_table_index_axis_list[start_index:]:
            next_plant_first_row = first_axis_table_index_axis_list.index(1, start_index)
        else:
            next_plant_first_row = len(first_axis_table_index_axis_list)
        current_plant_axes = first_axis_table_index_axis_list[current_plant_index:next_plant_first_row]
        index_plt_list.extend([plant_id for current_plant_axis in current_plant_axes])
        current_plant_index = next_plant_first_row
    return index_plt_list


def _create_index_axis_list(first_axis_table_plant_ids, cohort_probabilities):
    '''
    Create index_axis column.
    :Parameters:
        - `first_axis_table_plant_ids` : the plant indexes.
        - `cohort_probabilities` : the cohort probabilities.
    :Types:
        - `first_axis_table_plant_ids` : list
        - `cohort_probabilities` : dict
        
    :return: The index_axis column.
    :rtype: list
    '''
    index_axis_list = []
    for plant_id in first_axis_table_plant_ids:
        cohort_numbers = _find_child_cohort_numbers(cohort_probabilities)
        cohort_numbers.sort()
        index_axis_list.extend(cohort_numbers)
    return index_axis_list


def _find_child_cohort_numbers(cohort_probabilities, parent_cohort_number=-1, first_child_delay=2):
    '''
    Find (recursively) the child cohort numbers of a parent cohort, according to the cohort probabilities and 
    the parent cohort number. The main stem always exists.
    :Parameters:
        - `cohort_probabilities` : the cohort probabilities.
        - `parent_cohort_number` : the parent cohort number.
    :Types:
        - `cohort_probabilities` : dict
        - `parent_cohort_number` : int
        
    :return: The axes column.
    :rtype: list
    '''
    child_cohort_numbers = []
    first_possible_cohort_number = parent_cohort_number + first_child_delay
    if first_possible_cohort_number == 1:
        # The main stem always exists, then add it.
        child_cohort_numbers.append(first_possible_cohort_number)
        child_cohort_numbers.extend(_find_child_cohort_numbers(cohort_probabilities, 
                                                               first_possible_cohort_number))
    else:
        # Find the secondary stem children.
        for cohort_number_str, cohort_probability in cohort_probabilities.iteritems():
            cohort_number = int(cohort_number_str)
            if cohort_number >= first_possible_cohort_number:
                if cohort_probability >= random.random():
                    child_cohort_numbers.append(cohort_number)
                    child_cohort_numbers.extend(_find_child_cohort_numbers(cohort_probabilities, 
                                                                           cohort_number))
    return child_cohort_numbers
              
              
def _create_N_phyt_list(first_axis_table_index_axis_list, 
                        main_stem_leaves_number_probability_distribution,
                        secondary_stem_leaves_number_coefficients):
    '''
    Create the nff column.
    :Parameters:
        - `first_axis_table_index_axis_list` : the axes column.
        - `main_stem_leaves_number_probability_distribution` : the probability distribution of 
        the main stem leaves number.
        - `secondary_stem_leaves_number_coefficients` : the coefficients of the secondary stem leaves number.
    :Types:
        - `first_axis_table_index_axis_list` : list
        - `main_stem_leaves_number_probability_distribution` : dict
        - `secondary_stem_leaves_number_coefficients` : dict
        
    :return: The nff column.
    :rtype: list
    '''
    N_phyt_list = []
    main_stem_leaves_number = 0.0
    # for each plant...
    for cohort_number in first_axis_table_index_axis_list:
        # calculate the leaves number of each axis
        leaves_number_float = 0.0
        if cohort_number == 1:
            # It is the main stem, then the leaves number has to satisfy the probability distribution defined  
            # in main_stem_leaves_number_probability_distribution
            random_value = random.random()
            probabilities_sum = 0.0
            for leaves_number_str, leaves_probability in main_stem_leaves_number_probability_distribution.iteritems():
                probabilities_sum += leaves_probability
                if random_value <= probabilities_sum:
                    main_stem_leaves_number = float(leaves_number_str)
                    break
            leaves_number_float = main_stem_leaves_number
        else:
            # it is a secondary stem (i.e. a tiller)
            a_1 = secondary_stem_leaves_number_coefficients['a_1']
            a_2 = secondary_stem_leaves_number_coefficients['a_2']
            leaves_number_float = a_1* main_stem_leaves_number - a_2 * cohort_number
        fractional_part, integer_part = math.modf(leaves_number_float)
        if random.random() <= fractional_part:
            leaves_number_int = int(math.ceil(leaves_number_float))
        else:
            leaves_number_int = int(integer_part)
        N_phyt_list.append(leaves_number_int)
     
    return N_phyt_list


def _create_all_TT_phytomer1_list(first_axis_table_dataframe, emf_1_main_stem_standard_deviation, first_leaf_phen_table_dataframe):
    '''
    Create TT_em_phytomer1, TT_col_phytomer1, TT_sen_phytomer1 and TT_del_phytomer1.
    For each plant, define a delay of emergence, and for each axis add this delay to the first leaf development schedule.
    :Parameters:
        - `first_axis_table_dataframe` : the first axis table.
        - `emf_1_main_stem_standard_deviation` : the standard deviation used to calculate main stem emf_1 value.
        - `first_leaf_phen_table_dataframe` : the first leaf phen table dataframe.
    :Types:
        - `first_axis_table_dataframe` : pandas.DataFrame
        - `emf_1_main_stem_standard_deviation` : float
        - `first_leaf_phen_table_dataframe` : pandas.DataFrame
        
    :return: The TT_em_phytomer1, TT_col_phytomer1, TT_sen_phytomer1 and TT_del_phytomer1 data.
    :rtype: tuple of pandas.Series
    '''    
    sigma = emf_1_main_stem_standard_deviation
    sigma_div_2 = sigma / 2.0
    TT_em_phytomer1_series = pandas.Series(index=first_axis_table_dataframe.index)
    TT_col_phytomer1_series = pandas.Series(index=first_axis_table_dataframe.index)
    TT_sen_phytomer1_series = pandas.Series(index=first_axis_table_dataframe.index)
    TT_del_phytomer1_series = pandas.Series(index=first_axis_table_dataframe.index)

    for id_plt, first_axis_table_grouped_by_id_plt_dataframe in first_axis_table_dataframe.groupby('id_plt'):
        normal_distribution = random.normalvariate(0.0, sigma)
        while abs(normal_distribution) > sigma_div_2:
            normal_distribution = random.normalvariate(0.0, sigma)
        for id_phen, first_axis_table_grouped_by_id_plt_and_id_phen_dataframe in first_axis_table_grouped_by_id_plt_dataframe.groupby('id_phen'):
            current_first_leaf_phen_table_dataframe_row = first_leaf_phen_table_dataframe[first_leaf_phen_table_dataframe['id_phen']==id_phen]
            first_valid_index = current_first_leaf_phen_table_dataframe_row.first_valid_index()
            TT_em_phytomer1_series[first_axis_table_grouped_by_id_plt_and_id_phen_dataframe.index] = normal_distribution + current_first_leaf_phen_table_dataframe_row['TT_em_phytomer'][first_valid_index]
            TT_col_phytomer1_series[first_axis_table_grouped_by_id_plt_and_id_phen_dataframe.index] = normal_distribution + current_first_leaf_phen_table_dataframe_row['TT_col_phytomer'][first_valid_index]
            TT_sen_phytomer1_series[first_axis_table_grouped_by_id_plt_and_id_phen_dataframe.index] = normal_distribution + current_first_leaf_phen_table_dataframe_row['TT_sen_phytomer'][first_valid_index]
            TT_del_phytomer1_series[first_axis_table_grouped_by_id_plt_and_id_phen_dataframe.index] = normal_distribution + current_first_leaf_phen_table_dataframe_row['TT_del_phytomer'][first_valid_index]
                
    return TT_em_phytomer1_series, TT_col_phytomer1_series, TT_sen_phytomer1_series, TT_del_phytomer1_series  


def _create_id_dim_list(first_axis_table_index_axis_list, first_axis_table_N_phyt_list):
    '''
    Create id_dim column.
    :Parameters:
        - `first_axis_table_index_axis_list` : the axes column.
        - `first_axis_table_N_phyt_list` : the nff column.
    :Types:
        - `first_axis_table_index_axis_list` : list
        - `first_axis_table_N_phyt_list` : list
        
    :return: The id_dim column.
    :rtype: list
    '''
    return _create_ids(first_axis_table_index_axis_list, first_axis_table_N_phyt_list)


def _create_id_phen_list(first_axis_table_index_axis_list, first_axis_table_N_phyt_list):
    '''
    Create id_phen column.
    :Parameters:
        - `first_axis_table_index_axis_list` : the axes column.
        - `first_axis_table_N_phyt_list` : the nff column.
    :Types:
        - `first_axis_table_index_axis_list` : list
        - `first_axis_table_N_phyt_list` : list
        
    :return: The id_phen column.
    :rtype: list
    '''
    return _create_ids(first_axis_table_index_axis_list, first_axis_table_N_phyt_list)


def _create_ids(first_axis_table_index_axis_list, first_axis_table_N_phyt_list):
    '''
    Create id column.
    :Parameters:
        - `first_axis_table_index_axis_list` : the axes column.
        - `first_axis_table_N_phyt_list` : the nff column.
    :Types:
        - `first_axis_table_index_axis_list` : list
        - `first_axis_table_N_phyt_list` : list
        
    :return: The id column.
    :rtype: list
    '''
    id_list = []
    for i in range(len(first_axis_table_index_axis_list)):
        id_list.append(''.join([str(first_axis_table_index_axis_list[i]), str(first_axis_table_N_phyt_list[i]).zfill(2)]))
    return id_list


def _create_id_ear_list(first_axis_table_index_plt_list):
    '''
    Create id_ear column.
    :Parameters:
        - `first_axis_table_index_plt_list` : the plant ids column.
    :Types:
        - `first_axis_table_index_plt_list` : list
        
    :return: The id_ear column.
    :rtype: list
    '''
    return ['1' for plant_id in first_axis_table_index_plt_list]
    
    
def _create_TT_stop_axis_list(max_axes_number, min_axes_number, first_axis_table_TT_em_phytomer1, bolting_date, flowering_date):
    '''
    Create TT_stop_axis column.
    :Parameters:
        - `max_axes_number` : The maximum number of existing axes. Must be positive or null, and greater than min_axes_number.
        - `min_axes_number` : The minimum number of existing axes. Must be positive or null, and lesser than max_axes_number.
        - `first_axis_table_TT_em_phytomer1` : The TT_em_phytomer1 column.
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
    :Types:
        - `max_axes_number` : int
        - `min_axes_number` : int
        - `first_axis_table_TT_em_phytomer1` : list
        - `bolting_date` : int
        - `flowering_date` : int
        
    :return: The TT_stop_axis column.
    :rtype: pandas.Series
    '''
    
    assert max_axes_number >= 0 and min_axes_number >=0 and bolting_date >= 0 and flowering_date >= 0
    assert bolting_date <= flowering_date
    assert min_axes_number <= max_axes_number
    
    polynomial_coefficient_array = np.polyfit([flowering_date, bolting_date], [min_axes_number, max_axes_number], 1)
                
    remaining_axes_number = max_axes_number
    T_em_leaf1_tuples = zip(first_axis_table_TT_em_phytomer1[:], range(len(first_axis_table_TT_em_phytomer1)))
    T_em_leaf1_tuples.sort()
    T_stop_axis_tuples = []
    for tt in range(bolting_date, flowering_date + 1):
        simulated_axes_number = int(np.polyval(polynomial_coefficient_array, tt))
        axes_to_delete_number = remaining_axes_number - simulated_axes_number
        while axes_to_delete_number > 0:
            max_emf_1, axis_row_number = T_em_leaf1_tuples.pop()
            T_stop_axis_tuples.append((axis_row_number, tt))
            axes_to_delete_number -= 1
            remaining_axes_number -= 1
        if remaining_axes_number == 0:
            break 
    T_stop_axis_tuples.sort()
    T_stop_axis_row_number_list = [T_stop_axis_tuple[0] for T_stop_axis_tuple in T_stop_axis_tuples]
    TT_stop_axis_list = [T_stop_axis_tuple[1] for T_stop_axis_tuple in T_stop_axis_tuples]
    for i in range(len(first_axis_table_TT_em_phytomer1)):
        if i not in T_stop_axis_row_number_list:
            TT_stop_axis_list.insert(i, np.nan)
    return TT_stop_axis_list 
    
    
def _create_TT_del_axis_list(TT_stop_axis_series, delais_TT_stop_del_axis):
    '''
    Create TT_del_axis column.
    :Parameters:
        - `TT_stop_axis_series` : The TT_stop_axis column.
        - `delais_TT_stop_del_axis` : ???

    :Types:
        - `TT_stop_axis_series` : pandas.Series
        - `delais_TT_stop_del_axis` : int
        
    :return: The end column.
    :rtype: list
    '''
    return TT_stop_axis_series + delais_TT_stop_del_axis


def create_tillering_dynamic_dataframe(initial_date, bolting_date, flowering_date, plant_number, axis_table_dataframe, final_axes_number):
    '''
    Create the tillering dynamic dataframe. 
    :Parameters:
        - `initial_date` : the absolute phen table dataframe.
        - `bolting_date` : the bolting date
        - `flowering_date` : the flowering date
        - `plant_number` : the number of plants
        - `axis_table_dataframe` : the first axis table dataframe
        - `final_axes_number` : the final number of axes.
    :Types:
        - `initial_date` : int
        - `bolting_date` : int
        - `flowering_date` : int
        - `plant_number` : int
        - `axis_table_dataframe` : pandas.DataFrame
        - `final_axes_number` : int
        
    :return: the tillering dynamic dataframe.
    :rtype: pandas.Dataframe
    '''
    return pandas.DataFrame({'TT': [initial_date, bolting_date, flowering_date], 'NbrAxes': [plant_number, axis_table_dataframe.index.size, final_axes_number]}, columns=['TT', 'NbrAxes'])


