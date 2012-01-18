'''
This module provides functions to calculate AxisTable values.

Created on 28 nov. 2011

@author: cchambon
'''

import pandas
import random
import numpy
import math
try:
    from scipy import stats
except:
    from scipy import stats


def create_index_plt_list(plant_ids, index_axis_list):
    '''
    Create plant indexes column.
    :Parameters:
        - `plant_ids` : the plant indexes.
        - `index_axis_list` : the axes column.
    :Types:
        - `plant_ids` : list
        - `index_axis_list` : list
        
    :return: The plant indexes column.
    :rtype: list
    '''
    index_plt_list = []
    current_plant_index = 0
    for plant_id in plant_ids:
        start_index = current_plant_index + 1
        if 1 in index_axis_list[start_index:]:
            next_plant_first_row = index_axis_list.index(1, start_index)
        else:
            next_plant_first_row = len(index_axis_list)
        current_plant_axes = index_axis_list[current_plant_index:next_plant_first_row]
        index_plt_list.extend([plant_id for current_plant_axis in current_plant_axes])
        current_plant_index = next_plant_first_row
    return index_plt_list


def create_index_axis_list(plant_ids, 
                           cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.967, '6': 0.817, 
                                                 '7': 0.083, '8': 0.0, '9': 0.0, '10': 0.0}):
    '''
    Create axes column.
    :Parameters:
        - `plant_ids` : the plant indexes.
        - `cohort_probabilities` : the cohort probabilities.
    :Types:
        - `plant_ids` : list
        - `cohort_probabilities` : dict
        
    :return: The axes column.
    :rtype: list
    '''
    index_axis_list = []
    for plant_id in plant_ids:
        cohort_numbers = _find_child_cohort_numbers(cohort_probabilities)
        cohort_numbers.sort()
        index_axis_list.extend(cohort_numbers)
    return index_axis_list


def _find_child_cohort_numbers(cohort_probabilities, parent_cohort_number=-1):
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
    first_possible_cohort_number = parent_cohort_number + 2
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
              
              
def create_N_phyt_list(index_axis_list, 
                       main_stem_leaves_number_probability_distribution={'10': 0.145, '11': 0.818, '12': 0.036, 
                                                                         '13': 0.0, '14': 0.0},
                       secondary_stem_leaves_number_coefficients={'a_1': 0.9423, 'a_2': 0.555}):
    '''
    Create nff column.
    :Parameters:
        - `index_axis_list` : the axes column.
        - `main_stem_leaves_number_probability_distribution` : the probability distribution of 
        the main stem leaves number.
        - `secondary_stem_leaves_number_coefficients` : the coefficients of the secondary stem leaves number.
    :Types:
        - `index_axis_list` : list
        - `main_stem_leaves_number_probability_distribution` : dict
        - `secondary_stem_leaves_number_coefficients` : dict
        
    :return: The nff column.
    :rtype: list
    '''
    N_phyt_list = []
    main_stem_leaves_number = 0.0
    # for each plant...
    for cohort_number in index_axis_list:
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


def create_T_em_leaf1_list(index_axis_list, emf_1_main_stem_standard_deviation=30.0):
    '''
    Create emf_1 column.
    :Parameters:
        - `index_axis_list` : the axes column.
        - `emf_1_main_stem_standard_deviation` : the standard deviation used to calculate main stem emf_1 value.
    :Types:
        - `index_axis_list` : list
        - `emf_1_main_stem_standard_deviation` : float
        
    :return: The emf_1 column.
    :rtype: list
    '''
    mu=0.0
    sigma=emf_1_main_stem_standard_deviation
    T_em_leaf1_list = []
    for cohort_number in index_axis_list:
        emf_1 = 0.0
        if cohort_number == 1:
            # then it is the main stem
            emf_1 = random.normalvariate(mu, sigma)
            while abs(emf_1) > sigma / 2.0:
                emf_1 = random.normalvariate(mu, sigma)
        else:
            # it is a secondary stem
            emf_1 = None # TODO: will be modified.
        T_em_leaf1_list.append(emf_1)
    return T_em_leaf1_list


def create_id_dim_list(index_axis_list, N_phyt_list):
    '''
    Create id_dim column.
    :Parameters:
        - `index_axis_list` : the axes column.
        - `N_phyt_list` : the nff column.
    :Types:
        - `index_axis_list` : list
        - `N_phyt_list` : list
        
    :return: The id_dim column.
    :rtype: list
    '''
    return _create_ids(index_axis_list, N_phyt_list)


def create_id_phen_list(index_axis_list, N_phyt_list):
    '''
    Create id_phen column.
    :Parameters:
        - `index_axis_list` : the axes column.
        - `N_phyt_list` : the nff column.
    :Types:
        - `index_axis_list` : list
        - `N_phyt_list` : list
        
    :return: The id_phen column.
    :rtype: list
    '''
    return _create_ids(index_axis_list, N_phyt_list)


def _create_ids(index_axis_list, N_phyt_list):
    '''
    Create id column.
    :Parameters:
        - `index_axis` : the axes column.
        - `N_phyt_list` : the nff column.
    :Types:
        - `index_axis` : list
        - `N_phyt_list` : list
        
    :return: The id column.
    :rtype: list
    '''
    id_list = []
    for i in range(len(index_axis_list)):
        id_list.append(''.join([str(index_axis_list[i]), str(N_phyt_list[i]).zfill(2)]))
    return id_list


def create_id_ear_list(index_plt_list):
    '''
    Create id_ear column.
    :Parameters:
        - `index_plt_list` : the plant ids column.
    :Types:
        - `index_plt_list` : list
        
    :return: The id_ear column.
    :rtype: list
    '''
    return ['1' for plant_id in index_plt_list]
    
    
def create_T_stop_axis_list(max_axes_number, min_axes_number, T_em_leaf1_list, bolting_date=500, flowering_date=1000):
    '''
    Create end column.
    :Parameters:
        - `max_axes_number` : The maximum number of existing axes. Must be positive or null, and greater than min_axes_number.
        - `min_axes_number` : The minimum number of existing axes. Must be positive or null, and lesser than max_axes_number.
        - `T_em_leaf1_list` : The emf_1 column.
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
    :Types:
        - `max_axes_number` : int
        - `min_axes_number` : int
        - `T_em_leaf1_list` : list
        - `bolting_date` : int
        - `flowering_date` : int
        
    :return: The end column.
    :rtype: list
    '''
    
    assert max_axes_number >= 0 and min_axes_number >=0 and bolting_date >= 0 and flowering_date >= 0
    assert bolting_date < flowering_date
    assert min_axes_number < max_axes_number
    
    slope, intercept, r_value, p_value, std_err = stats.linregress([bolting_date, flowering_date],[max_axes_number, min_axes_number])              
    remaining_axes_number = max_axes_number
    T_em_leaf1_tuples = zip(T_em_leaf1_list[:], range(len(T_em_leaf1_list)))
    T_em_leaf1_tuples.sort()
    T_stop_axis_tuples = []
    for tt in range(bolting_date, flowering_date + 1):
        simulated_axes_number = int(slope * tt + intercept)
        axes_to_delete_number = remaining_axes_number - simulated_axes_number
        while axes_to_delete_number >= 0:
            max_emf_1, axis_row_number = T_em_leaf1_tuples.pop()
            T_stop_axis_tuples.append((axis_row_number, tt))
            axes_to_delete_number -= 1
            remaining_axes_number -= 1
        if remaining_axes_number == 0:
            break 
    T_stop_axis_tuples.sort()
    end_row_number = [T_stop_axis_tuple[0] for T_stop_axis_tuple in T_stop_axis_tuples]
    T_stop_axis_list = [T_stop_axis_tuple[1] for T_stop_axis_tuple in T_stop_axis_tuples]
    for i in range(len(T_em_leaf1_list)):
        if i not in end_row_number:
            T_stop_axis_list.insert(i, None) 
    return T_stop_axis_list 


if __name__ == "__main__":
    
    plant_ids = range(1,101)
    
    index_axis_list = create_index_axis_list(plant_ids)
    index_plt_list = create_index_plt_list(plant_ids, index_axis_list)
    N_phyt_list = create_N_phyt_list(index_axis_list)
    T_em_leaf1_list = create_T_em_leaf1_list(index_axis_list)

	# Remarque: avant de remplir la colonne TT_stop_axis il faut que la colonne TT_em_leaf1 soit totalement remplie (MB et Talles)

    T_stop_axis_list = create_T_stop_axis_list(len(index_axis_list), int(len(index_axis_list)/2), T_em_leaf1_list)
    id_dim_list = create_id_dim_list(index_axis_list, N_phyt_list)
    id_phen_list = create_id_phen_list(index_axis_list, N_phyt_list)
    id_ear_list = create_id_ear_list(index_plt_list)
    
    axis_array = numpy.array([index_plt_list, index_axis_list, N_phyt_list, T_stop_axis_list, id_dim_list, id_phen_list, id_ear_list, T_em_leaf1_list]).transpose()
    axis_table_dataframe = pandas.DataFrame(axis_array, columns=['id_plt', 'id_axis', 'N_phytomer', 'TT_stop_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_em_phytomer1'])
        
    import tempfile
    from openalea.core.path import path
    fitting_results_directory = path(tempfile.mkdtemp(suffix='_fitting_results'))
    axis_table_filepath = fitting_results_directory/'axis_table_random.csv'
    
    axis_table_dataframe.to_csv(axis_table_filepath, na_rep='NA', index=False)  
    
    print 'The results has been saved in %s' % fitting_results_directory
    
    import os
    os.system("start %s" % axis_table_filepath.abspath())

