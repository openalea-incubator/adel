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


def create_plant_ids_column(plant_ids, axes_column):
    '''
    Create plant indexes column.
    :Parameters:
        - `plant_ids` : the plant indexes.
        - `axes_column` : the axes column.
    :Types:
        - `plant_ids` : list
        - `axes_column` : list
        
    :return: The plant indexes column.
    :rtype: list
    '''
    plant_ids_column = []
    current_plant_index = 0
    for plant_id in plant_ids:
        start_index = current_plant_index + 1
        if 1 in axes_column[start_index:]:
            next_plant_first_row = axes_column.index(1, start_index)
        else:
            next_plant_first_row = len(axes_column)
        current_plant_axes = axes_column[current_plant_index:next_plant_first_row]
        plant_ids_column.extend([plant_id for current_plant_axis in current_plant_axes])
        current_plant_index = next_plant_first_row
    return plant_ids_column


def create_axes_column(plant_ids, 
                       cohort_probabilities={'3': 1.0, '4': 1.0, '5': 0.6, '6': 0.0, 
                                             '7': 0.0, '8': 0.0, '9': 0.0, '10': 0.0}):
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
    axes_column = []
    for plant_id in plant_ids:
        cohort_numbers = _find_child_cohort_numbers(cohort_probabilities)
        cohort_numbers.sort()
        axes_column.extend(cohort_numbers)
    return axes_column


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
              
              
def create_nff_column(axes_column, 
                      main_stem_leaves_number_probability_distribution={'10': 0.1, '11': 0.1, '12': 0.1, 
                                                                       '13': 0.1, '14': 0.6},
                      secondary_stem_leaves_number_coefficients={'a_1': 0.9423, 'a_2': 0.555}):
    '''
    Create nff column.
    :Parameters:
        - `axes_column` : the axes column.
        - `main_stem_leaves_number_probability_distribution` : the probability distribution of 
        the main stem leaves number.
        - `secondary_stem_leaves_number_coefficients` : the coefficients of the secondary stem leaves number.
    :Types:
        - `axes_column` : list
        - `main_stem_leaves_number_probability_distribution` : dict
        - `secondary_stem_leaves_number_coefficients` : dict
        
    :return: The nff column.
    :rtype: list
    '''
    nff_column = []
    main_stem_leaves_number = 0.0
    # for each plant...
    for cohort_number in axes_column:
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
        nff_column.append(leaves_number_int)
     
    return nff_column


def create_emf_1_column(axes_column, emf_1_main_stem_standard_deviation=30.0):
    '''
    Create emf_1 column.
    :Parameters:
        - `axes_column` : the axes column.
        - `emf_1_main_stem_standard_deviation` : the standard deviation used to calculate main stem emf_1 value.
    :Types:
        - `axes_column` : list
        - `emf_1_main_stem_standard_deviation` : float
        
    :return: The emf_1 column.
    :rtype: list
    '''
    mu=0.0
    sigma=emf_1_main_stem_standard_deviation
    emf_1_column = []
    for cohort_number in axes_column:
        emf_1 = 0.0
        if cohort_number == 1:
            # then it is the main stem
            emf_1 = random.normalvariate(mu, sigma)
            while abs(emf_1) > sigma / 2.0:
                emf_1 = random.normalvariate(mu, sigma)
        else:
            # it is a secondary stem
            emf_1 = None # TODO: will be modified.
        emf_1_column.append(emf_1)
    return emf_1_column


def create_id_dim_column(axes_column, nff_column):
    '''
    Create id_dim column.
    :Parameters:
        - `axes_column` : the axes column.
        - `nff_column` : the nff column.
    :Types:
        - `axes_column` : list
        - `nff_column` : list
        
    :return: The id_dim column.
    :rtype: list
    '''
    return _create_id_column(axes_column, nff_column)


def create_id_phen_column(axes_column, nff_column):
    '''
    Create id_phen column.
    :Parameters:
        - `axes_column` : the axes column.
        - `nff_column` : the nff column.
    :Types:
        - `axes_column` : list
        - `nff_column` : list
        
    :return: The id_phen column.
    :rtype: list
    '''
    return _create_id_column(axes_column, nff_column)


def _create_id_column(axes_column, nff_column):
    '''
    Create id column.
    :Parameters:
        - `axes_column` : the axes column.
        - `nff_column` : the nff column.
    :Types:
        - `axes_column` : list
        - `nff_column` : list
        
    :return: The id column.
    :rtype: list
    '''
    id_column = []
    for i in range(len(axes_column)):
        id_column.append(''.join([str(axes_column[i]), str(nff_column[i]).zfill(2)]))
    return id_column


def create_id_ear_column(plant_ids_column):
    '''
    Create id_ear column.
    :Parameters:
        - `plant_ids_column` : the plant ids column.
    :Types:
        - `plant_ids_column` : list
        
    :return: The id_ear column.
    :rtype: list
    '''
    return ['1' for plant_id in plant_ids_column]
    
    
def create_end_column(max_axes_number, min_axes_number, emf_1_column, bolting_date=500, flowering_date=1000):
    '''
    Create end column.
    :Parameters:
        - `max_axes_number` : The maximum number of existing axes. Must be positive or null, and greater than min_axes_number.
        - `min_axes_number` : The minimum number of existing axes. Must be positive or null, and lesser than max_axes_number.
        - `emf_1_column` : The emf_1 column.
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
    :Types:
        - `max_axes_number` : int
        - `min_axes_number` : int
        - `emf_1_column` : list
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
    emf_1_column_tuples = zip(emf_1_column[:], range(len(emf_1_column)))
    emf_1_column_tuples.sort()
    end_column_tuples = []
    for tt in range(bolting_date, flowering_date + 1):
        simulated_axes_number = int(slope * tt + intercept)
        axes_to_delete_number = remaining_axes_number - simulated_axes_number
        while axes_to_delete_number >= 0:
            max_emf_1, axis_row_number = emf_1_column_tuples.pop()
            end_column_tuples.append((axis_row_number, tt))
            axes_to_delete_number -= 1
            remaining_axes_number -= 1
        if remaining_axes_number == 0:
            break 
    end_column_tuples.sort()
    end_row_number = [end_column_tuple[0] for end_column_tuple in end_column_tuples]
    end_column = [end_column_tuple[1] for end_column_tuple in end_column_tuples]
    for i in range(len(emf_1_column)):
        if i not in end_row_number:
            end_column.insert(i, None) 
    return end_column 


if __name__ == "__main__":
    
    plant_ids = range(1,101)

    axes_column = create_axes_column(plant_ids)
    plant_ids_column = create_plant_ids_column(plant_ids, axes_column)
    nff_column = create_nff_column(axes_column)
    emf_1_column = create_emf_1_column(axes_column)
    end_column = create_end_column(len(axes_column), int(len(axes_column)/2), emf_1_column)
    id_dim_column = create_id_dim_column(axes_column, nff_column)
    id_phen_column = create_id_phen_column(axes_column, nff_column)
    id_ear_column = create_id_ear_column(plant_ids_column)
    
    axis_array = numpy.array([plant_ids_column, axes_column, nff_column, end_column, id_dim_column, id_phen_column, id_ear_column, emf_1_column]).transpose()
    
    import tempfile
    from openalea.core.path import path
    axis_table_directory = path(tempfile.mkdtemp(suffix='_axis_table_filling_test'))
    
    axis_dataframe = pandas.DataFrame(axis_array, columns=['plant', 'axe', 'nff', 'end', 'id_dim', 'id_phen', 'id_ear', 'emf_1'])
    axis_dataframe.to_csv(axis_table_directory/'axis_table_filling_test.csv', na_rep='NA', index=False)  
    
    print 'The results has been saved in %s' % axis_table_directory

