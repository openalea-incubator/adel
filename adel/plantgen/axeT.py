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
    Compute the following columns: *id_cohort_axis*, *id_plt*, *N_phytomer* and *id_phen*. 
           
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

    '''
    
    plant_ids = range(1,plant_number + 1)
    id_cohort_axis_list, id_axis_list = _gen_id_axis_list(plant_ids, decide_child_cohort_probabilities)
    id_plt_list = _gen_id_plt_list(plant_ids, id_cohort_axis_list)
    N_phytomer_list = _gen_N_phytomer_list(id_cohort_axis_list, 
                                           MS_leaves_number_probabilities, 
                                           params.SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS)
    HS_final_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_stop_axis_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_del_axis_list = [np.nan for i in range(len(id_cohort_axis_list))]
    id_dim_list = [np.nan for i in range(len(id_cohort_axis_list))]
    id_phen_list = _gen_id_phen_list(id_cohort_axis_list, N_phytomer_list)
    id_ear_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_em_phytomer1_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_col_phytomer1_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_sen_phytomer1_list = [np.nan for i in range(len(id_cohort_axis_list))]
    TT_del_phytomer1_list = [np.nan for i in range(len(id_cohort_axis_list))]
    axeT_array = np.array([id_plt_list, id_cohort_axis_list, N_phytomer_list, HS_final_list, TT_stop_axis_list, TT_del_axis_list, id_dim_list, id_phen_list, id_ear_list, TT_em_phytomer1_list, TT_col_phytomer1_list, TT_sen_phytomer1_list, TT_del_phytomer1_list]).transpose()
    axeT_array_df = pandas.DataFrame(axeT_array, columns=['id_plt', 'id_cohort_axis', 'N_phytomer', 'HS_final', 'TT_stop_axis', 'TT_del_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_em_phytomer1', 'TT_col_phytomer1', 'TT_sen_phytomer1', 'TT_del_phytomer1'])
    axeT_array_df.insert(2, 'id_axis', id_axis_list)
    return axeT_array_df


def create_axeT(axeT_tmp_dataframe, phenT_first_dataframe, dynT_dataframe, TT_bolting, TT_flag_leaf_ligulation, delais_TT_stop_del_axis, final_axes_density):
    '''
    Create the :ref:`axeT <axeT>` dataframe filling the *axeT_tmp* dataframe.
    
    :Parameters:
    
        - `axeT_tmp_dataframe` (:class:`pandas.DataFrame`) - the *axeT_tmp* dataframe.
        - `phenT_first_dataframe` (:class:`pandas.DataFrame`) - the :ref:`phenT_first <phenT_first>` dataframe.  
        - `dynT_dataframe` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe.
        - `TT_bolting` (:class:`float`) - date in thermal time at which the bolting starts.
        - `TT_flag_leaf_ligulation` (:class:`float`) - the thermal time of the flag leaf ligulation. 
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
          
    '''
    
    axeT_dataframe = axeT_tmp_dataframe.copy()
    (axeT_dataframe['TT_em_phytomer1'], 
     axeT_dataframe['TT_col_phytomer1'], 
     axeT_dataframe['TT_sen_phytomer1'],
     axeT_dataframe['TT_del_phytomer1']) = _gen_all_TT_phytomer1_list(axeT_tmp_dataframe, params.EMF_1_MS_STANDARD_DEVIATION, phenT_first_dataframe)
    axeT_dataframe['TT_stop_axis'] = tools.decide_time_of_death(axeT_tmp_dataframe.index.size, final_axes_density, axeT_dataframe['TT_em_phytomer1'].tolist(), TT_bolting, TT_flag_leaf_ligulation)
    axeT_dataframe['TT_del_axis'] = _gen_TT_del_axis_list(axeT_dataframe['TT_stop_axis'], delais_TT_stop_del_axis)
    axeT_dataframe['HS_final'] = _gen_HS_final_list(axeT_dataframe, dynT_dataframe)
    axeT_dataframe['id_dim'] = _gen_id_dim_list(axeT_dataframe['id_cohort_axis'], axeT_dataframe['HS_final'], axeT_dataframe['N_phytomer'])
    axeT_dataframe['id_ear'] = _gen_id_ear_list(axeT_dataframe['TT_stop_axis'])
    
    return axeT_dataframe
    

def _gen_id_plt_list(plant_ids, id_cohort_axis_list):
    '''Generate the *id_plt* column.'''
    id_plt_list = []
    current_plant_index = 0
    for plant_id in plant_ids:
        start_index = current_plant_index + 1
        if 1 in id_cohort_axis_list[start_index:]:
            next_plant_first_row = id_cohort_axis_list.index(1, start_index)
        else:
            next_plant_first_row = len(id_cohort_axis_list)
        current_plant_axes = id_cohort_axis_list[current_plant_index:next_plant_first_row]
        id_plt_list.extend([plant_id for current_plant_axis in current_plant_axes])
        current_plant_index = next_plant_first_row
    return id_plt_list


def _gen_id_axis_list(plant_ids, decide_child_cohort_probabilities):
    '''Generate the columns *id_axis* and *id_cohort_axis* .'''
    all_child_cohorts = []
    for plant_id in plant_ids:
        child_cohorts = tools.decide_child_cohorts(decide_child_cohort_probabilities, first_child_delay=params.FIRST_CHILD_DELAY)
        child_cohorts.sort()
        all_child_cohorts.extend(child_cohorts)
    all_child_cohorts_array = np.array(all_child_cohorts)
    cohort_numbers = all_child_cohorts_array[:, 0].astype(int).tolist()
    cohort_positions = all_child_cohorts_array[:, 1].tolist()
    return (cohort_numbers, cohort_positions)


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


def _gen_id_dim_list(id_cohort_axis_series, HS_final_series, N_phytomer_series):
    '''Generate the *id_dim* column.'''
    is_ear = pandas.Series(1, index=HS_final_series.index)
    is_ear[HS_final_series.dropna().index] = 0
    HS_final_series.fillna(N_phytomer_series, inplace=True)
    ceiled_HS_final_series = np.ceil(HS_final_series)
    zfilled_array = np.core.defchararray.zfill(np.char.mod('%d', ceiled_HS_final_series), 2)
    id_cohort_axis_str_array = np.char.mod('%d', id_cohort_axis_series)
    id_dim_array = np.core.defchararray.add(id_cohort_axis_str_array, zfilled_array)
    id_dim_array = np.core.defchararray.add(id_dim_array, np.char.mod('%d', is_ear)).astype(int)
    return id_dim_array.tolist()


def _gen_id_phen_list(id_cohort_axis_list, N_phyt_list):
    '''Generate the *id_phen* column.'''
    id_phen_list = []
    for i in range(len(id_cohort_axis_list)):
        id_phen_list.append(int(''.join([str(id_cohort_axis_list[i]), str(N_phyt_list[i]).zfill(2)])))
    return id_phen_list


def _gen_id_ear_list(TT_stop_axis):
    '''Generate the *id_ear* column.'''
    TT_stop_axis_series = pandas.Series(TT_stop_axis)
    id_ear = pandas.Series(1.0, index=TT_stop_axis_series.index)
    id_ear[TT_stop_axis_series.dropna().index] = np.nan
    return id_ear.tolist()
    
    
def _gen_TT_del_axis_list(TT_stop_axis_series, delais_TT_stop_del_axis):
    '''Generate the *TT_del_axis* column.'''
    return TT_stop_axis_series + delais_TT_stop_del_axis


def _gen_HS_final_list(axeT_dataframe, dynT_dataframe):
    '''Generate the *HS_final* column.'''
    HS_final_series = pandas.Series(index=axeT_dataframe.index)
    dynT_grouped = dynT_dataframe.groupby(['N_cohort', 'Nff'])
    for axeT_key, axeT_group in axeT_dataframe.groupby(['id_cohort_axis', 'N_phytomer']):
        dynT_group = dynT_grouped.get_group(axeT_key)
        current_a_cohort = dynT_group['a_cohort'][dynT_group.first_valid_index()]
        current_TT_col_0 = dynT_group['TT_col_0'][dynT_group.first_valid_index()]
        HS_final_series[axeT_group.index] = current_a_cohort * (axeT_group['TT_stop_axis'][axeT_group.index] - current_TT_col_0)
    index_to_modify = HS_final_series[HS_final_series > axeT_dataframe['N_phytomer']].index
    HS_final_series[index_to_modify] = axeT_dataframe['N_phytomer'][index_to_modify]
    return HS_final_series.tolist()


def create_tilleringT(initial_date, TT_bolting, TT_flag_leaf_ligulation, plant_number, axeT_tmp_dataframe, final_axes_density):
    '''
    Create the :ref:`tilleringT <tilleringT>` dataframe.
    
    :Parameters:
    
        - `initial_date` (:class:`int`) - the initial date.
        - `TT_bolting` (:class:`float`) - date in thermal time at which the bolting starts.
        - `TT_flag_leaf_ligulation` (:class:`float`) - the thermal time of the flag leaf ligulation.
        - `plant_number` (:class:`int`) - the number of plants. 
        - `axeT_tmp_dataframe` (:class:`pandas.Dataframe`) - the *axeT_tmp* dataframe.
        - `final_axes_density` (:class:`int`) - the final number of axes which have an ear, per square meter.
          
    :Returns:
        the :ref:`tilleringT <tilleringT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`

    '''
    
    return pandas.DataFrame({'TT': [initial_date, TT_bolting, TT_flag_leaf_ligulation], 'NbrAxes': [plant_number, axeT_tmp_dataframe.index.size, final_axes_density]}, columns=['TT', 'NbrAxes'])


def create_cohortT(theoretical_cohorts_cardinalities, id_cohort_axis):
    '''
    Create the :ref:`cohortT <cohortT>` dataframe.
    
    :Parameters:
    
        - `theoretical_cohorts_cardinalities` (:class:`dict`) - the theoretical 
          cardinalities of the cohorts. 
        - `id_cohort_axis` (:class:`pandas.Series`) - the *id_cohort_axis* column of 
          :ref:`axeT <axeT>`
          
    :Returns:
        the :ref:`cohortT <cohortT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
    
    '''
    theoretical_cohorts = theoretical_cohorts_cardinalities.keys()
    cohortT_dataframe = pandas.DataFrame(index=range(len(theoretical_cohorts)), 
                                         columns=['cohort', 
                                                  'theoretical_cardinality', 
                                                  'simulated_cardinality'])
    
    simulated_cardinalities = id_cohort_axis.astype(int).value_counts()
    idx = 0
    for theoretical_cohort, theoretical_cardinality in theoretical_cohorts_cardinalities.iteritems():
        cohortT_dataframe['cohort'][idx] = theoretical_cohort
        cohortT_dataframe['theoretical_cardinality'][idx] = theoretical_cardinality
        if theoretical_cohort in simulated_cardinalities:
            cohortT_dataframe['simulated_cardinality'][idx] = simulated_cardinalities[theoretical_cohort]
        else:
            cohortT_dataframe['simulated_cardinality'][idx] = 0
        idx += 1 
    cohortT_dataframe = cohortT_dataframe.sort_index(by='cohort')
    cohortT_dataframe.index = range(len(cohortT_dataframe))
    return cohortT_dataframe

