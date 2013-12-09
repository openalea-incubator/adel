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


def create_axeT_tmp(plants_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities):
    '''
    Create the *axeT_tmp* dataframe. 
    Compute the following columns: *id_plt*, *id_cohort*, *id_axis*, *N_phytomer_potential* and *id_phen*. 
           
    :Parameters:
    
        - `plants_number` (:class:`int`) - the number of plants. 
        - `decide_child_cohort_probabilities` (:class:`dict`) - the probabilities of the 
          child cohorts. 
        - `MS_leaves_number_probabilities` (:class:`dict`) - the probability 
          distribution of the main stem leaves number. 
          
    :Returns:
        the *axeT_tmp* dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`

    '''
    
    plant_ids = range(1, plants_number + 1)
    id_cohort_list, id_axis_list = _gen_id_axis_list(plant_ids, decide_child_cohort_probabilities)
    id_plt_list = _gen_id_plt_list(plant_ids, id_cohort_list)
    N_phytomer_potential_list = _gen_N_phytomer_potential_list(id_cohort_list, 
                                           MS_leaves_number_probabilities, 
                                           params.SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS)
    id_phen_list = _gen_id_phen_list(id_cohort_list, N_phytomer_potential_list)
    
    axeT_tmp = pandas.DataFrame(index=range(len(id_plt_list)),
                                   columns=['id_plt', 'id_cohort', 'id_axis', 'N_phytomer_potential', 'N_phytomer', 'HS_final', 'TT_stop_axis', 'TT_del_axis', 'id_dim', 'id_phen', 'id_ear', 'TT_app_phytomer1', 'TT_col_phytomer1', 'TT_sen_phytomer1', 'TT_del_phytomer1'],
                                   dtype=float)
    axeT_tmp['id_plt'] = id_plt_list
    axeT_tmp['id_cohort'] = id_cohort_list
    axeT_tmp['id_axis'] = id_axis_list
    axeT_tmp['N_phytomer_potential'] = N_phytomer_potential_list
    axeT_tmp['id_phen'] = id_phen_list

    return axeT_tmp


def create_axeT(axeT_tmp, phenT_first, dynT_, TT_regression_start, TT_flag_leaf_ligulation, delais_TT_stop_del_axis, number_of_ears):
    '''
    Create the :ref:`axeT <axeT>` dataframe filling the *axeT_tmp* dataframe.
    
    :Parameters:
    
        - `axeT_tmp` (:class:`pandas.DataFrame`) - the *axeT_tmp* dataframe.
        - `phenT_first` (:class:`pandas.DataFrame`) - the :ref:`phenT_first <phenT_first>` dataframe.  
        - `dynT_` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe.
        - `TT_regression_start` (:class:`float`) - thermal time at which the regression starts.
        - `TT_flag_leaf_ligulation` (:class:`float`) - the thermal time of the flag leaf ligulation. 
        - `delais_TT_stop_del_axis` (:class:`int`) - This variable represents the time in 
          thermal time between an axis stop growing and its disappearance (it 
          concerns only the axes that do not regress and which do not produce any 
          cob).
        - `number_of_ears` (:class:`int`) - the number of ears. 
          
    :Returns:
        the :ref:`axeT <axeT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
          
    '''
    
    axeT_ = axeT_tmp.copy()
    (axeT_['TT_app_phytomer1'], 
     axeT_['TT_col_phytomer1'], 
     axeT_['TT_sen_phytomer1'],
     axeT_['TT_del_phytomer1']) = _gen_all_TT_phytomer1_list(axeT_tmp, params.EMF_1_MS_STANDARD_DEVIATION, phenT_first)
    axeT_['TT_stop_axis'] = tools.decide_time_of_death(axeT_tmp.index.size, number_of_ears, axeT_['TT_app_phytomer1'].tolist(), TT_regression_start, TT_flag_leaf_ligulation)
    axeT_['id_ear'] = _gen_id_ear_list(axeT_['TT_stop_axis'])
    axeT_['TT_del_axis'] = _gen_TT_del_axis_list(axeT_['TT_stop_axis'], delais_TT_stop_del_axis)
    HS_final_series = _gen_HS_final_series(axeT_, dynT_)
    axeT_ = _remove_axes_without_leaf(axeT_, HS_final_series.index)
    axeT_['HS_final'] = HS_final_series.values
    axeT_['N_phytomer'] = _gen_N_phytomer(axeT_['HS_final'])
    axeT_['id_dim'] = _gen_id_dim_list(axeT_['id_cohort'], axeT_['N_phytomer'], axeT_['id_ear'])
    
    return axeT_
    

def _gen_id_plt_list(plant_ids, id_cohort_list):
    '''Generate the *id_plt* column.'''
    id_plt_list = []
    current_plant_index = 0
    for plant_id in plant_ids:
        start_index = current_plant_index + 1
        if 1 in id_cohort_list[start_index:]:
            next_plant_first_row = id_cohort_list.index(1, start_index)
        else:
            next_plant_first_row = len(id_cohort_list)
        current_plant_axes = id_cohort_list[current_plant_index:next_plant_first_row]
        id_plt_list.extend([plant_id for current_plant_axis in current_plant_axes])
        current_plant_index = next_plant_first_row
    return id_plt_list


def _gen_id_axis_list(plant_ids, decide_child_cohort_probabilities):
    '''Generate the columns *id_axis* and *id_cohort* .'''
    all_child_cohorts = []
    for plant_id in plant_ids:
        child_cohorts = tools.decide_child_cohorts(decide_child_cohort_probabilities, params.FIRST_CHILD_DELAY)
        child_cohorts.sort()
        all_child_cohorts.extend(child_cohorts)
    all_child_cohorts_array = np.array(all_child_cohorts)
    cohort_numbers = all_child_cohorts_array[:, 0].astype(int).tolist()
    cohort_positions = all_child_cohorts_array[:, 1].tolist()
    return (cohort_numbers, cohort_positions)


def _gen_N_phytomer_potential_list(id_cohort_list, 
                         MS_leaves_number_probabilities, 
                         secondary_stem_leaves_number_coefficients):
    '''Generate the *N_phytomer_potential* column.'''
    N_phytomer_potential_list = []
    MS_final_leaves_number = 0.0
    # for each plant...
    for cohort_number in id_cohort_list:
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
        N_phytomer_potential_list.append(leaves_number_int)
     
    return N_phytomer_potential_list


def _gen_N_phytomer(HS_final_series):
    '''Generate the *N_phytomer* column.'''
    return np.ceil(HS_final_series).astype(int)
    

def _gen_all_TT_phytomer1_list(axeT_tmp, emf_1_MS_standard_deviation, phenT_first):
    '''Generate the *TT_app_phytomer1*, *TT_col_phytomer1*, *TT_sen_phytomer1* and *TT_del_phytomer1* columns.
    For each plant, define a delay of appearance, and for each axis add this delay to the first leaf development schedule.'''
    sigma = emf_1_MS_standard_deviation
    sigma_div_2 = sigma / 2.0
    TT_app_phytomer1_series = pandas.Series(index=axeT_tmp.index)
    TT_col_phytomer1_series = pandas.Series(index=axeT_tmp.index)
    TT_sen_phytomer1_series = pandas.Series(index=axeT_tmp.index)
    TT_del_phytomer1_series = pandas.Series(index=axeT_tmp.index)

    for id_plt, axeT_tmp_grouped_by_id_plt in axeT_tmp.groupby('id_plt'):
        normal_distribution = random.normalvariate(0.0, sigma)
        while abs(normal_distribution) > sigma_div_2:
            normal_distribution = random.normalvariate(0.0, sigma)
        for id_phen, axeT_tmp_grouped_by_id_plt_and_id_phen in axeT_tmp_grouped_by_id_plt.groupby('id_phen'):
            current_row = phenT_first[phenT_first['id_phen']==id_phen]
            first_valid_index = current_row.first_valid_index()
            TT_app_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen.index] = normal_distribution + current_row['TT_app_phytomer'][first_valid_index]
            TT_col_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen.index] = normal_distribution + current_row['TT_col_phytomer'][first_valid_index]
            TT_sen_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen.index] = normal_distribution + current_row['TT_sen_phytomer'][first_valid_index]
            TT_del_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen.index] = normal_distribution + current_row['TT_del_phytomer'][first_valid_index]
                
    return TT_app_phytomer1_series, TT_col_phytomer1_series, TT_sen_phytomer1_series, TT_del_phytomer1_series  


def _gen_id_dim_list(id_cohort_series, N_phytomer_series, id_ear_series):
    '''Generate the *id_dim* column.'''
    is_ear = pandas.Series(0, index=id_ear_series.index)
    is_ear[id_ear_series.dropna().index] = 1
    zfilled_array = np.core.defchararray.zfill(np.char.mod('%d', N_phytomer_series), 2)
    id_cohort_str_array = np.char.mod('%d', id_cohort_series)
    id_dim_array = np.core.defchararray.add(id_cohort_str_array, zfilled_array)
    id_dim_array = np.core.defchararray.add(id_dim_array, np.char.mod('%d', is_ear)).astype(int)
    return id_dim_array.tolist()


def _gen_id_phen_list(id_cohort_list, N_phytomer_potential_list):
    '''Generate the *id_phen* column.'''
    id_phen_list = []
    for i in range(len(id_cohort_list)):
        id_phen_list.append(int(''.join([str(id_cohort_list[i]), str(N_phytomer_potential_list[i]).zfill(2)])))
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


def _gen_HS_final_series(axeT_, dynT_):
    '''Generate the *HS_final* column.'''
    HS_final_series = pandas.Series(index=axeT_.index)
    dynT_grouped = dynT_.groupby(['id_axis', 'N_phytomer_potential'])
    for axeT_key, axeT_group in axeT_.groupby(['id_axis', 'N_phytomer_potential']):
        dynT_group = dynT_grouped.get_group(axeT_key)
        current_a_cohort = dynT_group['a_cohort'][dynT_group.first_valid_index()]
        current_TT_col_0 = dynT_group['TT_col_0'][dynT_group.first_valid_index()]
        HS_final_series[axeT_group.index] = current_a_cohort * (axeT_group['TT_stop_axis'][axeT_group.index] - current_TT_col_0)
    index_to_modify = HS_final_series[HS_final_series > axeT_['N_phytomer_potential']].index
    HS_final_series[index_to_modify] = axeT_['N_phytomer_potential'][index_to_modify]
    HS_final_series.fillna(axeT_['N_phytomer_potential'], inplace=True)
    HS_final_series = HS_final_series.clip_lower(0.0)
    HS_final_series = HS_final_series[HS_final_series != 0.0]
    return HS_final_series


def _remove_axes_without_leaf(axeT_, index_to_keep):
    '''Remove the axes which do not have any leaf.'''
    axeT_ = axeT_.ix[index_to_keep]
    axeT_.index = range(len(axeT_))
    return axeT_
    

def create_tilleringT(TT_start, TT_regression_start, TT_flag_leaf_ligulation, plants_number, plants_density, number_of_axes, ears_density):
    '''
    Create the :ref:`tilleringT <tilleringT>` dataframe.
    
    :Parameters:
    
        - `TT_start` (:class:`int`) - the thermal time at which the growth starts.
        - `TT_regression_start` (:class:`float`) - the thermal time at which the regression starts.
        - `TT_flag_leaf_ligulation` (:class:`float`) - the thermal time of the flag leaf ligulation.
        - `plants_number` (:class:`int`) - the number of plants to simulate.
        - `plants_density` (:class:`int`) - the number of plants that are present 
          after loss due to bad emergence, early death..., per square meter.
        - `number_of_axes` (:class:`int`) - the number of simulated axes.
        - `ears_density` (:class:`int`) - the number of ears per square meter.
          
    :Returns:
        the :ref:`tilleringT <tilleringT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`

    '''
    axes_density = number_of_axes / float(plants_number) * plants_density 
    return pandas.DataFrame({'TT': [TT_start, TT_regression_start, TT_flag_leaf_ligulation], 'axes_density': [plants_density, axes_density, ears_density]}, columns=['TT', 'axes_density'])


def create_cardinalityT(theoretical_cohort_cardinalities, theoretical_axis_cardinalities, simulated_cohorts_axes):
    '''
    Create the :ref:`cardinalityT <cardinalityT>` dataframe.
    
    :Parameters:
    
        - `theoretical_cohort_cardinalities` (:class:`dict`) - the theoretical 
          cardinalities of the cohorts. 
        - `theoretical_axis_cardinalities` (:class:`dict`) - the theoretical 
          cardinalities of the axes. 
        - `simulated_cohorts_axes` (:class:`pandas.DataFrame`) - the *id_cohort* 
          and *id_axis* columns of :ref:`axeT <axeT>`.
          
    :Returns:
        the :ref:`cardinalityT <cardinalityT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
    
    '''
    simulated_cohort_cardinalities = simulated_cohorts_axes['id_cohort'].value_counts().to_dict()
    simulated_axis_cardinalities = simulated_cohorts_axes.groupby(['id_cohort', 'id_axis']).size().to_dict()
    cardinalityT = pandas.DataFrame(index=range(len(theoretical_axis_cardinalities)), 
                                    columns=['id_cohort', 
                                             'id_axis',
                                             'theoretical_cohort_cardinality', 
                                             'simulated_cohort_cardinality',
                                             'theoretical_axis_cardinality',
                                             'simulated_axis_cardinality'])
    idx = 0
    for (id_cohort, id_axis), theoretical_axis_cardinality in theoretical_axis_cardinalities.iteritems():
        cardinalityT['id_cohort'][idx] = id_cohort
        cardinalityT['id_axis'][idx] = id_axis
        cardinalityT['theoretical_cohort_cardinality'][idx] = theoretical_cohort_cardinalities[id_cohort]
        cardinalityT['theoretical_axis_cardinality'][idx] = theoretical_axis_cardinality
        if id_cohort in simulated_cohort_cardinalities:
            cardinalityT['simulated_cohort_cardinality'][idx] = simulated_cohort_cardinalities[id_cohort]
        else:
            cardinalityT['simulated_cohort_cardinality'][idx] = 0
        if (id_cohort, id_axis) in simulated_axis_cardinalities:
            cardinalityT['simulated_axis_cardinality'][idx] = simulated_axis_cardinalities[(id_cohort, id_axis)]
        else:
            cardinalityT['simulated_axis_cardinality'][idx] = 0
        idx += 1 
    cardinalityT[['theoretical_cohort_cardinality', 
                  'theoretical_axis_cardinality',
                  'simulated_cohort_cardinality',
                  'simulated_axis_cardinality']] = cardinalityT[['theoretical_cohort_cardinality', 
                                                                 'theoretical_axis_cardinality',
                                                                 'simulated_cohort_cardinality',
                                                                 'simulated_axis_cardinality']].astype(float)
    cardinalityT['id_cohort'] = cardinalityT['id_cohort'].astype(int)                                                              
    cardinalityT.sort(['id_cohort', 'id_axis'], inplace=True)
    cardinalityT.index = range(len(cardinalityT))
    return cardinalityT

