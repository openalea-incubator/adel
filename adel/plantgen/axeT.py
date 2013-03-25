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
'''
Provides functions to construct the *axeT_tmp*, the :ref:`axeT <axeT>` and the :ref:`tilleringT <tilleringT>` dataframes.

The :ref:`axeT <axeT>` and the :ref:`tilleringT <tilleringT>` dataframes are described in the User Guide 
(see :ref:`adel_user`).

Authors: Mariem Abichou, Camille Chambon, Bruno Andrieu
'''

import math

import random
import numpy as np
import pandas

from adel.plantgen import tools, params


def create_axeT_tmp(plant_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities):
    '''
    Create the *axeT_tmp* dataframe. 
    Compute the following columns: *id_cohort_axis*, *id_plt*, *N_phytomer*, *id_dim*, *id_phen*, *id_ear*. 
           
    :Parameters:
    
        - `plant_number` (:class:`int`) - the number of plants. 
        - `decide_child_cohort_probabilities` (:class:`dict`) - the probabilities of the 
          child cohorts. 
        - `MS_leaves_number_probabilities` (:class:`dict`) - the probability 
          distribution of the main stem leaves number. 
          
    :Returns:
        the *axeT_tmp* dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`

    .. warning:: the type of the arguments is checked as follows:

         .. list-table::
             :widths: 10 50
             :header-rows: 1
        
             * - Argument
               - Type
             * - *plant_number* 
               - :class:`int`
             * - *decide_child_cohort_probabilities* 
               - :class:`dict`
             * - *MS_leaves_number_probabilities* 
               - :class:`dict`

    '''
    assert isinstance(plant_number, int)
    assert isinstance(decide_child_cohort_probabilities, dict)
    assert isinstance(MS_leaves_number_probabilities, dict)
    
    plant_ids = range(1,plant_number + 1)
    id_cohort_axis_list = _gen_id_cohort_axis_list(plant_ids, decide_child_cohort_probabilities)
    index_plt_list = _gen_index_plt_list(plant_ids, id_cohort_axis_list)
    N_phytomer_list = _gen_N_phytomer_list(id_cohort_axis_list, 
                                       MS_leaves_number_probabilities, 
                                       params.secondary_stem_leaves_number_coefficients)
    TT_stop_axis_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_del_axis_list = [np.nan for i in range(len(id_cohort_axis_list))]
    id_dim_list = _gen_id_dim_list(id_cohort_axis_list, N_phytomer_list)
    id_phen_list = _gen_id_phen_list(id_cohort_axis_list, N_phytomer_list)
    id_ear_list = _gen_id_ear_list(index_plt_list)
    TT_em_phytomer1_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_col_phytomer1_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_sen_phytomer1_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_del_phytomer1_list = [np.nan for i in range(len(id_cohort_axis_list))]
    axeT_array = np.array([index_plt_list, id_cohort_axis_list, N_phytomer_list, TT_stop_axis_list, TT_del_axis_list, id_dim_list, id_phen_list, id_ear_list, TT_em_phytomer1_list, TT_col_phytomer1_list, TT_sen_phytomer1_list, TT_del_phytomer1_list]).transpose()
    return pandas.DataFrame(axeT_array, columns=['id_plt', 'id_cohort_axis', 'N_phytomer', 'TT_stop_axis', 'TT_del_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_em_phytomer1', 'TT_col_phytomer1', 'TT_sen_phytomer1', 'TT_del_phytomer1'], dtype=float)


def create_axeT(axeT_tmp_dataframe, phenT_first_dataframe, TT_bolting, TT_flowering, delais_TT_stop_del_axis, final_axes_density):
    '''
    Create the :ref:`axeT <axeT>` dataframe filling the *axeT_tmp* dataframe.
    
    :Parameters:
    
        - `axeT_tmp_dataframe` (:class:`pandas.DataFrame`) - the *axeT_tmp* dataframe.
        - `phenT_first_dataframe` (:class:`pandas.DataFrame`) - the :ref:`phenT_first <phenT_first>` dataframe.  
        - `TT_bolting` (:class:`float`) - date in thermal time at which the bolting starts.
        - `TT_flowering` (:class:`float`) - the flowering date. 
        - `delais_TT_stop_del_axis` (:class:`int`) - This variable represents the time in 
          thermal time between an axis stop growing and its disappearance (it 
          concerns only the axes that do not regress and which do not produce any 
          cob).
        - `final_axes_density` (:class:`int`) - the final number of axes which have an ear, 
          per square meter.
          
    :Returns:
        the :ref:`axeT <axeT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`

    .. warning:: the type of the arguments is checked as follows:

         .. list-table::
             :widths: 10 50
             :header-rows: 1
        
             * - Argument
               - Type
             * - *axeT_tmp_dataframe* 
               - :class:`pandas.DataFrame`
             * - *phenT_first_dataframe* 
               - :class:`pandas.DataFrame`
             * - *TT_bolting* 
               - :class:`float`
             * - *TT_flowering* 
               - :class:`float`
             * - *delais_TT_stop_del_axis* 
               - :class:`int`
             * - *final_axes_density* 
               - :class:`int`
          
    '''
    assert isinstance(axeT_tmp_dataframe, pandas.DataFrame)
    assert isinstance(phenT_first_dataframe, pandas.DataFrame)
    assert isinstance(TT_bolting, float)
    assert isinstance(TT_flowering, float)
    assert isinstance(delais_TT_stop_del_axis, int)
    assert isinstance(final_axes_density, int)
    axeT_dataframe = axeT_tmp_dataframe.copy()
    (axeT_dataframe['TT_em_phytomer1'], 
     axeT_dataframe['TT_col_phytomer1'], 
     axeT_dataframe['TT_sen_phytomer1'],
     axeT_dataframe['TT_del_phytomer1']) = _gen_all_TT_phytomer1_list(axeT_tmp_dataframe, params.emf_1_MS_standard_deviation, phenT_first_dataframe)
    axeT_dataframe['TT_stop_axis'] = tools.decide_time_of_death(axeT_tmp_dataframe.index.size, final_axes_density, axeT_dataframe['TT_em_phytomer1'].tolist(), TT_bolting, TT_flowering)
    axeT_dataframe['TT_del_axis'] = _gen_TT_del_axis_list(axeT_dataframe['TT_stop_axis'], delais_TT_stop_del_axis)
    
    return axeT_dataframe
    

def _gen_index_plt_list(plant_ids, id_cohort_axis_list):
    '''Generate the *id_plt* column.'''
    index_plt_list = []
    current_plant_index = 0
    for plant_id in plant_ids:
        start_index = current_plant_index + 1
        if 1 in id_cohort_axis_list[start_index:]:
            next_plant_first_row = id_cohort_axis_list.index(1, start_index)
        else:
            next_plant_first_row = len(id_cohort_axis_list)
        current_plant_axes = id_cohort_axis_list[current_plant_index:next_plant_first_row]
        index_plt_list.extend([plant_id for current_plant_axis in current_plant_axes])
        current_plant_index = next_plant_first_row
    return index_plt_list


def _gen_id_cohort_axis_list(plant_ids, decide_child_cohort_probabilities):
    '''Generate the *id_cohort_axis* column.'''
    id_cohort_axis_list = []
    for plant_id in plant_ids:
        cohort_numbers = tools.decide_child_cohorts(decide_child_cohort_probabilities, first_child_delay=params.first_child_delay)
        cohort_numbers.sort()
        id_cohort_axis_list.extend(cohort_numbers)
    return id_cohort_axis_list


def _gen_N_phytomer_list(id_cohort_axis_list, 
                         MS_leaves_number_probabilities, 
                         secondary_stem_leaves_number_coefficients):
    '''Generate the *N_phytomer* column.'''
    N_phytomer_list = []
    MS_final_leaves_number = 0.0
    # for each plant...
    for cohort_number in id_cohort_axis_list:
        # calculate the leaves number of each axis
        leaves_number_float = 0.0
        if cohort_number == 1:
            # It is the main stem, then the leaves number has to satisfy the probability distribution defined  
            # in MS_leaves_number_probabilities
            MS_final_leaves_number = tools.calculate_MS_final_leaves_number(MS_leaves_number_probabilities)
            leaves_number_float = MS_final_leaves_number
        else:
            # it is a secondary stem (i.e. a tiller)
            leaves_number_float = tools.calculate_tiller_final_leaves_number(MS_final_leaves_number, cohort_number, secondary_stem_leaves_number_coefficients)
        fractional_part, integer_part = math.modf(leaves_number_float)
        if random.random() <= fractional_part:
            leaves_number_int = int(math.ceil(leaves_number_float))
        else:
            leaves_number_int = int(integer_part)
        N_phytomer_list.append(leaves_number_int)
     
    return N_phytomer_list


def _gen_all_TT_phytomer1_list(axeT_tmp_dataframe, emf_1_MS_standard_deviation, phenT_first_dataframe):
    '''Generate the *TT_em_phytomer1*, *TT_col_phytomer1*, *TT_sen_phytomer1* and *TT_del_phytomer1* columns.
    For each plant, define a delay of emergence, and for each axis add this delay to the first leaf development schedule.'''
    sigma = emf_1_MS_standard_deviation
    sigma_div_2 = sigma / 2.0
    TT_em_phytomer1_series = pandas.Series(index=axeT_tmp_dataframe.index)
    TT_col_phytomer1_series = pandas.Series(index=axeT_tmp_dataframe.index)
    TT_sen_phytomer1_series = pandas.Series(index=axeT_tmp_dataframe.index)
    TT_del_phytomer1_series = pandas.Series(index=axeT_tmp_dataframe.index)

    for id_plt, axeT_tmp_grouped_by_id_plt_dataframe in axeT_tmp_dataframe.groupby('id_plt'):
        normal_distribution = random.normalvariate(0.0, sigma)
        while abs(normal_distribution) > sigma_div_2:
            normal_distribution = random.normalvariate(0.0, sigma)
        for id_phen, axeT_tmp_grouped_by_id_plt_and_id_phen_dataframe in axeT_tmp_grouped_by_id_plt_dataframe.groupby('id_phen'):
            current_row = phenT_first_dataframe[phenT_first_dataframe['id_phen']==id_phen]
            first_valid_index = current_row.first_valid_index()
            TT_em_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen_dataframe.index] = normal_distribution + current_row['TT_em_phytomer'][first_valid_index]
            TT_col_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen_dataframe.index] = normal_distribution + current_row['TT_col_phytomer'][first_valid_index]
            TT_sen_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen_dataframe.index] = normal_distribution + current_row['TT_sen_phytomer'][first_valid_index]
            TT_del_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen_dataframe.index] = normal_distribution + current_row['TT_del_phytomer'][first_valid_index]
                
    return TT_em_phytomer1_series, TT_col_phytomer1_series, TT_sen_phytomer1_series, TT_del_phytomer1_series  


def _gen_id_dim_list(id_cohort_axis_list, N_phyt_list):
    '''Generate the *id_dim* column.'''
    return _gen_ids(id_cohort_axis_list, N_phyt_list)


def _gen_id_phen_list(id_cohort_axis_list, N_phyt_list):
    '''Generate the *id_phen* column.'''
    return _gen_ids(id_cohort_axis_list, N_phyt_list)


def _gen_ids(id_cohort_axis_list, N_phyt_list):
    '''Generate an *id_\** column.'''
    id_list = []
    for i in range(len(id_cohort_axis_list)):
        id_list.append(''.join([str(id_cohort_axis_list[i]), str(N_phyt_list[i]).zfill(2)]))
    return id_list


def _gen_id_ear_list(index_plt_list):
    '''Generate the *id_ear* column.'''
    return ['1' for plant_id in index_plt_list]
    
    
def _gen_TT_del_axis_list(TT_stop_axis_series, delais_TT_stop_del_axis):
    '''Construct the *TT_del_axis* column.'''
    return TT_stop_axis_series + delais_TT_stop_del_axis


def create_tilleringT(initial_date, TT_bolting, TT_flowering, plant_number, axeT_tmp_dataframe, final_axes_density):
    '''
    Create the :ref:`tilleringT <tilleringT>` dataframe.
    
    :Parameters:
    
        - `initial_date` (:class:`int`) - the initial date.
        - `TT_bolting` (:class:`float`) - date in thermal time at which the bolting starts.
        - `TT_flowering` (:class:`float`) - the flowering date.
        - `plant_number` (:class:`int`) - the number of plants. 
        - `axeT_tmp_dataframe` (:class:`pandas.Dataframe`) - the *axeT_tmp* dataframe.
        - `final_axes_density` (:class:`int`) - the final number of axes which have an ear, per square meter.
          
    :Returns:
        the :ref:`tilleringT <tilleringT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`

    .. warning:: the type of the arguments is checked as follows:

         .. list-table::
             :widths: 10 50
             :header-rows: 1
        
             * - Argument
               - Type
             * - *initial_date* 
               - :class:`int`
             * - *TT_bolting* 
               - :class:`float`
             * - *TT_flowering* 
               - :class:`float`
             * - *plant_number* 
               - :class:`int`
             * - *axeT_tmp_dataframe* 
               - :class:`pandas.DataFrame`  
             * - *final_axes_density* 
               - :class:`int`
    
    '''
    assert isinstance(initial_date, int)
    assert isinstance(TT_bolting, float)
    assert isinstance(TT_flowering, float)
    assert isinstance(plant_number, int)
    assert isinstance(axeT_tmp_dataframe, pandas.DataFrame)
    assert isinstance(final_axes_density, int)
    return pandas.DataFrame({'TT': [initial_date, TT_bolting, TT_flowering], 'NbrAxes': [plant_number, axeT_tmp_dataframe.index.size, final_axes_density]}, columns=['TT', 'NbrAxes'])


def create_cohortT(plant_number, decide_child_cohort_probabilities, id_cohort_axes):
    '''
    Create the :ref:`cohortT <cohortT>` dataframe.
    
    :Parameters:
    
        - `plant_number` (:class:`int`) - the number of plants. 
        - `decide_child_cohort_probabilities` (:class:`dict`) - the probabilities of the child 
          cohorts. 
        - `id_cohort_axes` (:class:`pandas.Series`) - the *id_cohort_axis* column of 
          :ref:`axeT <axeT>`
          
    :Returns:
        the :ref:`cohortT <cohortT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`

    .. warning:: the type of the arguments is checked as follows:

         .. list-table::
             :widths: 10 50
             :header-rows: 1
        
             * - Argument
               - Type
             * - *plant_number* 
               - :class:`int`
             * - *decide_child_cohort_probabilities* 
               - :class:`dict`
             * - *id_cohort_axes* 
               - :class:`pandas.Series`
    
    '''
    assert isinstance(plant_number, int)
    assert isinstance(decide_child_cohort_probabilities, dict)
    assert isinstance(id_cohort_axes, pandas.Series)
    
    decide_cohort_probabilities = decide_child_cohort_probabilities.copy()
    # the cohort '1' always exists, so its probability is 1.0.
    decide_cohort_probabilities['1'] = 1.0
    possible_cohorts = np.array(decide_cohort_probabilities.keys()).astype(int)
    decide_cohort_probability_values = np.array(decide_cohort_probabilities.values())
    possible_child_cohorts = np.array(decide_child_cohort_probabilities.keys()).astype(int)
    decide_child_cohort_probability_values = np.array(decide_child_cohort_probabilities.values())
    possible_parent_cohorts = possible_child_cohorts - params.first_child_delay
    commons = np.where(np.intersect1d(possible_cohorts, possible_parent_cohorts))
    possible_parent_cohorts = possible_cohorts[commons]
    decide_parent_cohort_probability_values = decide_cohort_probability_values[commons]
    
    cohortT_dataframe = pandas.DataFrame(index=range(len(decide_cohort_probabilities)), 
                                         columns=['cohort', 
                                                  'theoretical_cardinality', 
                                                  'simulated_cardinality'])
    theoretical_probabilities_series = pandas.Series(index=cohortT_dataframe.index)
    simulated_cardinalities = id_cohort_axes.astype(int).value_counts()
    idx = 0
    for (cohort_str, decide_probability) in decide_cohort_probabilities.iteritems():
        cohort_int = int(cohort_str)
        cohortT_dataframe['cohort'][idx] = cohort_int
        if cohort_str == '1':
            theoretical_probabilities_series[idx] = decide_cohort_probabilities['1']
        else:
            first_possible_parent = cohort_int - params.first_child_delay
            curr_possible_parent_indexes = np.where(possible_parent_cohorts <= first_possible_parent)
            curr_decide_parent_cohort_probabilities = decide_parent_cohort_probability_values[curr_possible_parent_indexes]
            theoretical_probabilities_series[idx] = (curr_decide_parent_cohort_probabilities * decide_probability).sum()
        if cohort_int in simulated_cardinalities:
            cohortT_dataframe['simulated_cardinality'][idx] = simulated_cardinalities[cohort_int]
        else:
            cohortT_dataframe['simulated_cardinality'][idx] = 0
        idx += 1 
    cohortT_dataframe['theoretical_cardinality'] = theoretical_probabilities_series * plant_number
    cohortT_dataframe = cohortT_dataframe.sort_index(by='cohort')
    return cohortT_dataframe
        
