# -*- python -*-
#
#       Adel.Fit
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
import random
import math

import numpy as np

def find_child_cohort_numbers(cohort_probabilities, parent_cohort_number=-1, first_child_delay=2):
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
        child_cohort_numbers.extend(find_child_cohort_numbers(cohort_probabilities, 
                                                               first_possible_cohort_number))
    else:
        # Find the secondary stem children.
        for cohort_number_str, cohort_probability in cohort_probabilities.iteritems():
            cohort_number = int(cohort_number_str)
            if cohort_number >= first_possible_cohort_number:
                if cohort_probability >= random.random():
                    child_cohort_numbers.append(cohort_number)
                    child_cohort_numbers.extend(find_child_cohort_numbers(cohort_probabilities, 
                                                                           cohort_number))
    return child_cohort_numbers


def create_N_phyt_list(first_axis_table_index_axis_list, 
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
    MS_final_leaves_number = 0.0
    # for each plant...
    for cohort_number in first_axis_table_index_axis_list:
        # calculate the leaves number of each axis
        leaves_number_float = 0.0
        if cohort_number == 1:
            # It is the main stem, then the leaves number has to satisfy the probability distribution defined  
            # in main_stem_leaves_number_probability_distribution
            MS_final_leaves_number = _fit_MS_final_leaves_number(main_stem_leaves_number_probability_distribution)
            leaves_number_float = MS_final_leaves_number
        else:
            # it is a secondary stem (i.e. a tiller)
            leaves_number_float = _fit_tiller_final_leaves_number(MS_final_leaves_number, cohort_number, secondary_stem_leaves_number_coefficients)
        fractional_part, integer_part = math.modf(leaves_number_float)
        if random.random() <= fractional_part:
            leaves_number_int = int(math.ceil(leaves_number_float))
        else:
            leaves_number_int = int(integer_part)
        N_phyt_list.append(leaves_number_int)
     
    return N_phyt_list


def _fit_MS_final_leaves_number(main_stem_leaves_number_probability_distribution):
    random_value = random.random()
    probabilities_sum = 0.0
    for leaves_number_str, leaves_probability in main_stem_leaves_number_probability_distribution.iteritems():
        probabilities_sum += leaves_probability
        if random_value <= probabilities_sum:
            MS_final_leaves_number = float(leaves_number_str)
            break
    return MS_final_leaves_number


def _fit_tiller_final_leaves_number(MS_final_leaves_number, cohort_number, secondary_stem_leaves_number_coefficients):
    a_1 = secondary_stem_leaves_number_coefficients['a_1']
    a_2 = secondary_stem_leaves_number_coefficients['a_2']
    return a_1* MS_final_leaves_number - a_2 * cohort_number
    

def dead_or_alive_decision(max_axes_number, min_axes_number, first_axis_table_TT_em_phytomer1, bolting_date, flowering_date):
    '''
    Algorithme pour la decision de l'arret de croissance d'une talle en fonction d'une loi de decroissance d'effectif 
    sur la population globale.
    (previous name: Create TT_stop_axis column).
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

