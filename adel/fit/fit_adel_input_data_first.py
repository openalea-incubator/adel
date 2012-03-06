'''
This module is a facade for the first step of Adel input data fitting .

Created on Feb 7, 2012

@author: cchambon
'''  
import numpy
import pandas

import adel.fit.axis_table_fitting as axis_table_fitting
import adel.fit.dim_table_fitting as dim_table_fitting

def fit_adel_input_data_first(plant_number=100, 
                              cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.967, '6': 0.817, '7': 0.083, '8': 0.0, '9': 0.0, '10': 0.0}, 
                              main_stem_leaves_number_probability_distribution={'10': 0.145, '11': 0.818, '12': 0.036, '13': 0.0, '14': 0.0}, 
                              bolting_date=500, 
                              flowering_date=1000):
    '''
    Fit the axis table partially, initialize the parameters table and initialize the dim table.
    :Parameters:
        - `plant_number` : the number of plants.
        - `cohort_probabilities` : the cohort probabilities.
        - `main_stem_leaves_number_probability_distribution` : the probability distribution of the main stem leaves number.
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
    :Types:
        - `plant_number` : int.
        - `cohort_probabilities` : dict.
        - `main_stem_leaves_number_probability_distribution` : dict
        - `bolting_date` : int
        - `flowering_date` : int

    :return: The partially fitted axis table, the initialized dim table and the initialized parameters table.
    :rtype: a tuple of pandas.DataFrame
    '''    
    # Create and fit AxisTable
    axis_table_dataframe = axis_table_fitting.fit_axis_table_first(plant_number, cohort_probabilities, main_stem_leaves_number_probability_distribution, bolting_date, flowering_date)
    # Initialize the parameters table
    first_parameters_table_dataframe = fit_user_parameters_first(axis_table_dataframe['id_phen'].tolist())
    # Initialize DimTable
    first_dim_table_dataframe = dim_table_fitting.fit_dim_table_first(first_parameters_table_dataframe)

    return axis_table_dataframe, first_dim_table_dataframe, first_parameters_table_dataframe


def fit_user_parameters_first(first_axis_table_id_phen_list):
    '''
    Initialize the parameters table.
    :Parameters:
        - `first_axis_table_id_phen_list` : the number of plants.

    :Types:
        - `first_axis_table_id_phen_list` : int.

    :return: The initialized parameters table.
    :rtype: pandas.DataFrame
    ''' 
    id_phen_without_duplicate_list = list(set(first_axis_table_id_phen_list))
    N_cohort = [float(str(int(id_phen))[:-2]) for id_phen in id_phen_without_duplicate_list]
    axis_frequency_list = _create_axis_frequency_list(first_axis_table_id_phen_list, id_phen_without_duplicate_list)
    Nff = [float(str(int(id_phen))[-2:]) for id_phen in id_phen_without_duplicate_list]
    a_cohort_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_col_0_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_HS_break_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    TT_HS_NFF_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    dTT_MS_cohort_list = [numpy.nan for i in range(len(id_phen_without_duplicate_list))]
    parameters_table_array = numpy.array([N_cohort, id_phen_without_duplicate_list, axis_frequency_list, Nff, a_cohort_list, TT_col_0_list, TT_HS_break_list, TT_HS_NFF_list, dTT_MS_cohort_list]).transpose()
    unsorted_parameters_table_dataframe = pandas.DataFrame(parameters_table_array, columns=['N_cohort', 'id_axis', 'frequency', 'Nff', 'a_cohort', 'TT_col_0', 'TT_col_break', 'TT_col_nff', 'dTT_MS_cohort'], dtype=float)
    sorted_parameters_table_dataframe = pandas.DataFrame(columns=['N_cohort', 'id_axis', 'frequency', 'Nff', 'a_cohort', 'TT_col_0', 'TT_col_break', 'TT_col_nff', 'dTT_MS_cohort'], dtype=float)
    for name, group in unsorted_parameters_table_dataframe.groupby('N_cohort'):
        sorted_group = group.sort_index(by='frequency', ascending=False)
        sorted_parameters_table_dataframe = sorted_parameters_table_dataframe.append(sorted_group)
    return sorted_parameters_table_dataframe


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

