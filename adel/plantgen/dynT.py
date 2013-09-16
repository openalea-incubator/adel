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
Provides functions to construct the *dynT_tmp* and the :ref:`dynT <dynT>` dataframes.

The :ref:`dynT <dynT>` dataframe is described in the User Guide (see :ref:`adel_user`).

Authors: Mariem Abichou, Camille Chambon, Bruno Andrieu
'''

import numpy as np
import pandas

from adel.plantgen import params, tools


def create_dynT_tmp(axeT_tmp):
    '''
    Create the *dynT_tmp* dataframe.

    :Parameters:
        - `axeT_tmp` (:class:`pandas.DataFrame`) - the *axeT_tmp* dataframe. 

    :Returns:
        The *dynT_tmp* dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    '''
    groups = axeT_tmp.groupby(['id_axis', 'id_cohort', 'N_phytomer_potential']).groups
    keys_array = np.array(groups.keys())
    cardinalities = pandas.DataFrame(np.array(groups.values())).applymap(np.size)
    # initialize the values of the other columns to NaN
    dynT_tmp = pandas.DataFrame(index=range(len(groups)), 
                                columns=['id_axis', 'id_cohort', 'cardinality', 'N_phytomer_potential', 'a_cohort', 'TT_col_0', 'TT_col_break', 'TT_col_N_phytomer_potential', 'dTT_MS_cohort', 'n0', 'n1', 'n2', 't0', 't1', 'hs_t1', 'a', 'c', 'RMSE_gl'],
                                dtype=float)
    
    # set the columns 'id_axis', 'id_cohort', 'cardinality' and 'N_phytomer_potential'
    dynT_tmp['id_axis'] = keys_array[:, 0]
    dynT_tmp['id_cohort'] = keys_array[:, 1].astype(float).astype(int)
    dynT_tmp['cardinality'] = cardinalities
    dynT_tmp['N_phytomer_potential'] = keys_array[:, 2].astype(float).astype(int)
    
    # nested sort of dynT_tmp: first by 'id_axis' in ascending order, second 
    # by 'cardinality' in descending order.
    dynT_tmp.sort(['id_axis', 'cardinality'], ascending=[1, 0], inplace=True)
    # reinitialize the index
    dynT_tmp.index = range(dynT_tmp.index.size)
    return dynT_tmp


def create_dynT(dynT_tmp, 
                GL_number, 
                decimal_elongated_internode_number,
                leaf_number_delay_MS_cohort=params.LEAF_NUMBER_DELAY_MS_COHORT):
    '''
    Create the :ref:`dynT <dynT>` dataframe.
    
    :Parameters:
         - `dynT_tmp` (:class:`pandas.DataFrame`) - the *dynT_tmp* dataframe.
         - `GL_number` (:class:`dict` of :class:`float`::class:`float`) - the GL decimal number measured at 
           several thermal times (including the senescence end).
         - `decimal_elongated_internode_number` (:class:`float`) - the number of 
           elongated internodes.
         - `leaf_number_delay_MS_cohort` (:class:`dict` of :class:`int`::class:`float`) - the delays between 
           the emergence of the main stem and the emergence of each cohort.

    :Returns:
        The :ref:`dynT <dynT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * in *dynT_tmp*, the column *id_axis*, *id_phen*, *cardinality* 
          and *N_phytomer_potential* must be completely filled, i.e. they must not contain 
          any NA value.
    
    ''' 
    # in 'dynT_tmp', check that the columns 'id_axis', 'cardinality' and 
    # 'N_phytomer_potential' are non-NA.
    if dynT_tmp['id_axis'].count() != dynT_tmp['id_axis'].size:
        raise tools.InputError("dynT_tmp['id_axis'] contains NA values")
    if dynT_tmp['cardinality'].count() != dynT_tmp['cardinality'].size:
        raise tools.InputError("dynT_tmp['cardinality'] contains NA values")
    if dynT_tmp['N_phytomer_potential'].count() != dynT_tmp['N_phytomer_potential'].size:
        raise tools.InputError("dynT_tmp['N_phytomer_potential'] contains NA values")
    # get all main stem rows
    MS = dynT_tmp[dynT_tmp['id_axis'] == 'MS']
    
    # get the row of the most frequent main stem
    most_frequent_MS = MS.ix[0:0]
    # for this row, fill the columns referring to the dynamic of the green leaves
    most_frequent_MS = _gen_most_frequent_MS_GL_dynamic(most_frequent_MS, decimal_elongated_internode_number, GL_number)
    
    # get the rows of all main stems except the most frequent one
    other_MS = MS.ix[1:]
    # for these rows, fill the columns referring to the dynamic of Haun Stage
    other_MS = _gen_other_MS_HS_dynamic(most_frequent_MS, other_MS)
    # and fill the columns referring to the dynamic of the green leaves
    other_MS = _gen_other_MS_GL_dynamic(most_frequent_MS, other_MS)
    
    # get the rows corresponding to the tiller axes
    tiller_axes = dynT_tmp[dynT_tmp['id_axis'] != 'MS']
    # extract the rows corresponding to the most frequent tiller axes
    grouped = tiller_axes.groupby('id_axis')
    most_frequent_tiller_axes = []
    for id_axis, group_indexes in grouped.groups.iteritems():
        most_frequent_tiller_axes.append(tiller_axes.ix[group_indexes[0:1]])
    # concatenate these rows in one dataframe ; 'most_frequent_tiller_axes' is  
    # now a pandas.DataFrame (and is not a list anymore)
    most_frequent_tiller_axes = pandas.concat(most_frequent_tiller_axes)
    # in this new dataframe, fill the columns referring to the dynamic of Haun Stage
    most_frequent_tiller_axes = _gen_most_frequent_tiller_axes_HS_dynamic(most_frequent_MS, most_frequent_tiller_axes)
    # and fill the columns referring to the dynamic of the green leaves
    most_frequent_tiller_axes = _gen_most_frequent_tiller_axes_GL_dynamic(most_frequent_MS, most_frequent_tiller_axes)
    # extract the rows corresponding to all tiller axes except the most frequent ones
    other_tiller_axes = tiller_axes.drop(most_frequent_tiller_axes.index)
    other_tiller_axes = _gen_other_tiller_axes_HS_dynamic(most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes)
    # and fill the columns referring to the dynamic of the green leaves
    other_tiller_axes = _gen_other_tiller_axes_GL_dynamic(most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes)
    
    dynT_ = pandas.concat([most_frequent_MS, other_MS, most_frequent_tiller_axes, other_tiller_axes])
    
    dynT_.sort(['id_axis', 'cardinality'], ascending=[1, 0], inplace=True)
    # reinitialize the index
    dynT_.index = range(dynT_.index.size)
    
    return dynT_
    

def _gen_most_frequent_MS_GL_dynamic(most_frequent_MS, decimal_elongated_internode_number, GL_number):
    '''
    Create a copy of *most_frequent_MS*, fill this copy by calculating the 
    parameters which describe the dynamic of the green leaves, and return this 
    copy.  
    '''
    # calculation of t1
    most_frequent_MS = most_frequent_MS.copy()
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_MS['t1'] = most_frequent_MS['TT_col_0'] + decimal_elongated_internode_number / most_frequent_MS['a_cohort']
    else: # bilinear mode
        HS_break_0 = most_frequent_MS['a_cohort'][0] * (most_frequent_MS['TT_col_break'][0] - most_frequent_MS['TT_col_0'][0])
        a2_0 = (most_frequent_MS['N_phytomer_potential'][0] - HS_break_0) / (most_frequent_MS['TT_col_N_phytomer_potential'][0] - most_frequent_MS['TT_col_break'][0])
        if decimal_elongated_internode_number < HS_break_0:
            most_frequent_MS['t1'] = most_frequent_MS['TT_col_0'] + decimal_elongated_internode_number / most_frequent_MS['a_cohort']
        else:
            most_frequent_MS['t1'] = (decimal_elongated_internode_number - HS_break_0) / a2_0 + most_frequent_MS['TT_col_break']
    # calculation of hs_t1
    most_frequent_MS['hs_t1'] = most_frequent_MS['a_cohort'] * (most_frequent_MS['t1'] - most_frequent_MS['TT_col_0'])
    # calculation of t0
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_MS['t0'] = most_frequent_MS['TT_col_0'] + most_frequent_MS['n0'] / most_frequent_MS['a_cohort']
    else: # bilinear mode
        HS_break = most_frequent_MS['a_cohort'] * (most_frequent_MS['TT_col_break'] - most_frequent_MS['TT_col_0'])
        a2 = (most_frequent_MS['N_phytomer_potential'] - most_frequent_MS['HS_break']) / (most_frequent_MS['TT_col_N_phytomer_potential'] - most_frequent_MS['TT_col_break'])
        n0_smaller_than_HS_break_indexes = most_frequent_MS[most_frequent_MS['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = most_frequent_MS[most_frequent_MS['n0'] >= HS_break].index
        most_frequent_MS['t0'][n0_smaller_than_HS_break_indexes] = most_frequent_MS['TT_col_0'][n0_smaller_than_HS_break_indexes] + most_frequent_MS['n0'][n0_smaller_than_HS_break_indexes] / most_frequent_MS['a_cohort'][n0_smaller_than_HS_break_indexes]
        most_frequent_MS['t0'][n0_greater_than_HS_break_indexes] = (most_frequent_MS['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + most_frequent_MS['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # calculation of c 
    most_frequent_MS['c'] = -((most_frequent_MS['N_phytomer_potential'] - decimal_elongated_internode_number) - (most_frequent_MS['n2'] - most_frequent_MS['n1'])) / (most_frequent_MS['TT_col_N_phytomer_potential'] - most_frequent_MS['t1'])
    # calculation of a
    TT_col_N_phytomer_potential_0 = most_frequent_MS['TT_col_N_phytomer_potential'][0]
    n2_0 = most_frequent_MS['n2'][0]
    TT = np.array([TT_col_N_phytomer_potential_0] + GL_number.keys()) - TT_col_N_phytomer_potential_0
    GL = np.array([n2_0] + GL_number.values())
    fixed_coefs = [0.0, most_frequent_MS['c'][0], n2_0]
    a_starting_estimate = -4.0e-9
    most_frequent_MS['a'][0], most_frequent_MS['RMSE_gl'][0] = \
    tools.fit_poly(TT, GL, fixed_coefs, a_starting_estimate)
    return most_frequent_MS


def _gen_other_MS_HS_dynamic(most_frequent_MS, other_MS):
    '''
    Create a copy of *other_MS*, fill this copy by calculating the 
    parameters which describe the dynamic of the Haun Stage, and return this 
    copy.
    '''
    other_MS = other_MS.copy()  
    if other_MS['TT_col_0'].count() != other_MS['TT_col_0'].size:
        # calculation of TT_col_0
        other_MS['TT_col_0'] = most_frequent_MS['TT_col_0'][0]
        # calculation of TT_col_break
        other_MS['TT_col_break'] = most_frequent_MS['TT_col_break'][0]
        # calculation of dTT_MS_cohort
        other_MS['dTT_MS_cohort'] = most_frequent_MS['dTT_MS_cohort'][0] + (other_MS['id_cohort'] - most_frequent_MS['id_cohort'][0]) / (4 * most_frequent_MS['a_cohort'][0])
        # calculation of TT_col_N_phytomer_potential
        other_MS['TT_col_N_phytomer_potential'] = other_MS['dTT_MS_cohort'] + most_frequent_MS['TT_col_N_phytomer_potential'][0]
        # calculation of a_cohort
        if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
            other_MS['a_cohort'] = other_MS['N_phytomer_potential'] / (other_MS['TT_col_N_phytomer_potential'] - other_MS['TT_col_0'])
        else: # bilinear mode
            other_MS['a_cohort'] = most_frequent_MS['a_cohort'][0] 
    return other_MS


def _gen_other_MS_GL_dynamic(most_frequent_MS, other_MS):
    '''
    Create a copy of *other_MS*, fill this copy by calculating the 
    parameters which describe the dynamic of the green leaves, and return this 
    copy.  
    '''
    other_MS = other_MS.copy()
    # calculation of n1
    if other_MS['n1'].count() != other_MS['n1'].size:
        other_MS['n1'] = most_frequent_MS['n1'][0] * other_MS['N_phytomer_potential'] / most_frequent_MS['N_phytomer_potential'][0]
    # calculation of t1
    other_MS['t1'] = most_frequent_MS['t1'][0] + other_MS['dTT_MS_cohort']
    # calculation of hs_t1
    other_MS['hs_t1'] = other_MS['a_cohort'] * (other_MS['t1'] - other_MS['TT_col_0'])
    # calculation of n0 
    if other_MS['n0'].count() != other_MS['n0'].size:
        other_MS['n0'] = most_frequent_MS['n0'][0] * other_MS['N_phytomer_potential'] / most_frequent_MS['N_phytomer_potential'][0]
    # calculation of n2
    if other_MS['n2'].count() != other_MS['n2'].size:
        other_MS['n2'] = most_frequent_MS['n2'][0] * other_MS['N_phytomer_potential'] / most_frequent_MS['N_phytomer_potential'][0]
    # calculation of t0
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        other_MS['t0'] = other_MS['TT_col_0'] + other_MS['n0'] / other_MS['a_cohort']
    else: # bilinear mode
        HS_break = other_MS['a_cohort'] * (other_MS['TT_col_break'] - other_MS['TT_col_0'])
        a2 = (other_MS['N_phytomer_potential'] - other_MS['HS_break']) / (other_MS['TT_col_N_phytomer_potential'] - other_MS['TT_col_break'])
        n0_smaller_than_HS_break_indexes = other_MS[other_MS['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = other_MS[other_MS['n0'] >= HS_break].index
        other_MS['t0'][n0_smaller_than_HS_break_indexes] = other_MS['TT_col_0'][n0_smaller_than_HS_break_indexes] + other_MS['n0'][n0_smaller_than_HS_break_indexes] / other_MS['a_cohort'][n0_smaller_than_HS_break_indexes]
        other_MS['t0'][n0_greater_than_HS_break_indexes] = (other_MS['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + other_MS['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # calculation of c
    other_MS['c'] = most_frequent_MS['c'][0] * other_MS['n2'] / most_frequent_MS['n2'][0]
    # calculation of a
    other_MS['a'] = most_frequent_MS['a'][0] * other_MS['n2'] / most_frequent_MS['n2'][0]
    # calculation of RMSE_gl
    other_MS['RMSE_gl'] = most_frequent_MS['RMSE_gl'][0]
    return other_MS


def _gen_most_frequent_tiller_axes_HS_dynamic(most_frequent_MS, most_frequent_tiller_axes):
    '''
    Create a copy of *most_frequent_tiller_axes*, fill this copy by calculating the 
    parameters which describe the dynamic of the Haun Stage, and return this 
    copy.
    '''
    most_frequent_tiller_axes = most_frequent_tiller_axes.copy()
    # calculation of TT_col_break
    most_frequent_tiller_axes['TT_col_break'] = most_frequent_MS['TT_col_break'][0]
    without_nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes.dropna(subset=['TT_col_0']).index
    nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes.index - without_nan_most_frequent_tiller_axis_indexes
    # calculation of TT_col_0
    cohorts = most_frequent_tiller_axes['id_axis'][nan_most_frequent_tiller_axis_indexes].values
    leaf_number_delay_MS_cohorts = np.array([params.LEAF_NUMBER_DELAY_MS_COHORT[cohort] for cohort in cohorts if cohort in params.LEAF_NUMBER_DELAY_MS_COHORT])
    most_frequent_tiller_axes['TT_col_0'][nan_most_frequent_tiller_axis_indexes] = most_frequent_MS['TT_col_0'][0] + (leaf_number_delay_MS_cohorts / most_frequent_MS['a_cohort'][0])
    # dTT_MS_cohort is set by the user. Thus there is nothing to do.
    # calculation of TT_col_N_phytomer_potential
    most_frequent_tiller_axes['TT_col_N_phytomer_potential'][nan_most_frequent_tiller_axis_indexes] = most_frequent_tiller_axes['dTT_MS_cohort'][nan_most_frequent_tiller_axis_indexes] + most_frequent_MS['TT_col_N_phytomer_potential'][0]
    # calculation of a_cohort
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_tiller_axes['a_cohort'][nan_most_frequent_tiller_axis_indexes] = most_frequent_tiller_axes['N_phytomer_potential'][nan_most_frequent_tiller_axis_indexes] / (most_frequent_tiller_axes['TT_col_N_phytomer_potential'][nan_most_frequent_tiller_axis_indexes] - most_frequent_tiller_axes['TT_col_0'][nan_most_frequent_tiller_axis_indexes])
    else: # bilinear mode
        most_frequent_tiller_axes['a_cohort'][nan_most_frequent_tiller_axis_indexes] = most_frequent_MS['a_cohort'][0]
    return most_frequent_tiller_axes


def _gen_most_frequent_tiller_axes_GL_dynamic(most_frequent_MS, most_frequent_tiller_axes):
    '''
    Create a copy of *most_frequent_tiller_axes*, fill this copy by calculating the 
    parameters which describe the dynamic of the green leaves, and return this 
    copy.  
    '''
    most_frequent_tiller_axes = most_frequent_tiller_axes.copy()
    without_nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes.dropna(subset=['n1']).index
    nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes.index - without_nan_most_frequent_tiller_axis_indexes
    # calculation of n1
    most_frequent_tiller_axes['n1'][nan_most_frequent_tiller_axis_indexes] = most_frequent_MS['n1'][0]
    # calculation of t1
    most_frequent_tiller_axes['t1'] = most_frequent_MS['t1'][0] + most_frequent_tiller_axes['dTT_MS_cohort']
    # calculation of hs_t1
    most_frequent_tiller_axes['hs_t1'] = most_frequent_tiller_axes['a_cohort'] * (most_frequent_tiller_axes['t1'] - most_frequent_tiller_axes['TT_col_0'])
    # calculation of n0
    n1_hs_t1 = most_frequent_tiller_axes[['n1', 'hs_t1']].ix[nan_most_frequent_tiller_axis_indexes]
    if n1_hs_t1.index.size != 0:
        most_frequent_tiller_axes['n0'][nan_most_frequent_tiller_axis_indexes] = n1_hs_t1.apply(np.min, 1)                      
    # calculation of n2
    most_frequent_tiller_axes['n2'][nan_most_frequent_tiller_axis_indexes] = most_frequent_MS['n2'][0] * params.N2_MS_DIV_N2_COHORT
    # calculation of t0
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_tiller_axes['t0'] = most_frequent_tiller_axes['TT_col_0'] + most_frequent_tiller_axes['n0'] / most_frequent_tiller_axes['a_cohort']
    else: # bilinear mode
        HS_break = most_frequent_tiller_axes['a_cohort'] * (most_frequent_tiller_axes['TT_col_break'] - most_frequent_tiller_axes['TT_col_0'])
        a2 = (most_frequent_tiller_axes['N_phytomer_potential'] - most_frequent_tiller_axes['HS_break']) / (most_frequent_tiller_axes['TT_col_N_phytomer_potential'] - most_frequent_tiller_axes['TT_col_break'])
        n0_smaller_than_HS_break_indexes = most_frequent_tiller_axes[most_frequent_tiller_axes['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = most_frequent_tiller_axes[most_frequent_tiller_axes['n0'] >= HS_break].index
        most_frequent_tiller_axes['t0'][n0_smaller_than_HS_break_indexes] = most_frequent_tiller_axes['TT_col_0'][n0_smaller_than_HS_break_indexes] + most_frequent_tiller_axes['n0'][n0_smaller_than_HS_break_indexes] / most_frequent_tiller_axes['a_cohort'][n0_smaller_than_HS_break_indexes]
        most_frequent_tiller_axes['t0'][n0_greater_than_HS_break_indexes] = (most_frequent_tiller_axes['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + most_frequent_tiller_axes['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # calculation of c
    most_frequent_tiller_axes['c'] = most_frequent_MS['c'][0] * most_frequent_tiller_axes['n2'] / most_frequent_MS['n2'][0]
    # calculation of a
    most_frequent_tiller_axes['a'] = most_frequent_MS['a'][0] * most_frequent_tiller_axes['n2'] / most_frequent_MS['n2'][0]
    # calculation of RMSE_gl
    most_frequent_tiller_axes['RMSE_gl'] = most_frequent_MS['RMSE_gl'][0]
    return most_frequent_tiller_axes
    

def _gen_other_tiller_axes_HS_dynamic(most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes):
    '''
    Create a copy of *other_tiller_axes*, fill this copy by calculating the 
    parameters which describe the dynamic of the Haun Stage, and return this 
    copy.
    '''
    other_tiller_axes = other_tiller_axes.copy()
    # calculation of TT_col_break
    other_tiller_axes['TT_col_break'] = most_frequent_MS['TT_col_break'][0]    
    without_nan_other_tiller_axis_indexes = other_tiller_axes.dropna(subset=['TT_col_0']).index
    nan_other_tiller_axis_indexes = other_tiller_axes.index - without_nan_other_tiller_axis_indexes
    for name, group in other_tiller_axes.ix[nan_other_tiller_axis_indexes].groupby('id_axis'):
        most_frequent_tiller_axis_idx = most_frequent_tiller_axes[most_frequent_tiller_axes['id_axis'] == name].first_valid_index()
        # calculation of TT_col_0
        other_tiller_axes['TT_col_0'][group.index] = most_frequent_tiller_axes['TT_col_0'][most_frequent_tiller_axis_idx]
        # calculation of dTT_MS_cohort
        other_tiller_axes['dTT_MS_cohort'][group.index] = \
            most_frequent_tiller_axes['dTT_MS_cohort'][most_frequent_tiller_axis_idx] \
            + (group['id_cohort'] - most_frequent_tiller_axes['id_cohort'][most_frequent_tiller_axis_idx]) \
            / (4 * most_frequent_MS['a_cohort'][0])
        if most_frequent_MS['TT_col_break'][0] != 0.0: 
            # calculation of a_cohort in bilinear mode
            other_tiller_axes['a_cohort'][group.index] = most_frequent_tiller_axes['a_cohort'][most_frequent_tiller_axis_idx]
    # calculation of TT_col_N_phytomer_potential
    other_tiller_axes['TT_col_N_phytomer_potential'][nan_other_tiller_axis_indexes] = other_tiller_axes['dTT_MS_cohort'][nan_other_tiller_axis_indexes] + most_frequent_MS['TT_col_N_phytomer_potential'][0]
    # calculation of a_cohort in linear mode
    if most_frequent_MS['TT_col_break'][0] == 0.0:
        other_tiller_axes['a_cohort'][nan_other_tiller_axis_indexes] = other_tiller_axes['N_phytomer_potential'][nan_other_tiller_axis_indexes] / (other_tiller_axes['TT_col_N_phytomer_potential'][nan_other_tiller_axis_indexes] - other_tiller_axes['TT_col_0'][nan_other_tiller_axis_indexes])
    return other_tiller_axes
    

def _gen_other_tiller_axes_GL_dynamic(most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes):
    '''
    Create a copy of *other_tiller_axes*, fill this copy by calculating the 
    parameters which describe the dynamic of the green leaves, and return this 
    copy.  
    '''
    other_tiller_axes = other_tiller_axes.copy()
    without_nan_other_tiller_axis_indexes = other_tiller_axes.dropna(subset=['n1']).index
    nan_other_tiller_axis_indexes = other_tiller_axes.index - without_nan_other_tiller_axis_indexes
    for name, group in other_tiller_axes.ix[nan_other_tiller_axis_indexes].groupby('id_axis'):
        most_frequent_tiller_axis_idx = most_frequent_tiller_axes[most_frequent_tiller_axes['id_axis'] == name].first_valid_index()
        # calculation ofn1
        other_tiller_axes['n1'][group.index] = most_frequent_tiller_axes['n1'][most_frequent_tiller_axis_idx]
        # calculation of n2
        other_tiller_axes['n2'][group.index] = most_frequent_tiller_axes['n2'][most_frequent_tiller_axis_idx]
    # calculation of t1
    other_tiller_axes['t1'] = most_frequent_MS['t1'][0] + other_tiller_axes['dTT_MS_cohort']
    # calculation of hs_t1
    other_tiller_axes['hs_t1'] = other_tiller_axes['a_cohort'] * (other_tiller_axes['t1'] - other_tiller_axes['TT_col_0'])
    # calculation of n0
    n1_hs_t1 = other_tiller_axes[['n1', 'hs_t1']].ix[nan_other_tiller_axis_indexes]
    if n1_hs_t1.index.size != 0:
        other_tiller_axes['n0'][nan_other_tiller_axis_indexes] = n1_hs_t1.apply(np.min, 1)                      
    # calculation of t0
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        other_tiller_axes['t0'] = other_tiller_axes['TT_col_0'] + other_tiller_axes['n0'] / other_tiller_axes['a_cohort']
    else: # bilinear mode
        HS_break = other_tiller_axes['a_cohort'] * (other_tiller_axes['TT_col_break'] - other_tiller_axes['TT_col_0'])
        a2 = (other_tiller_axes['N_phytomer_potential'] - other_tiller_axes['HS_break']) / (other_tiller_axes['TT_col_N_phytomer_potential'] - other_tiller_axes['TT_col_break'])
        n0_smaller_than_HS_break_indexes = other_tiller_axes[other_tiller_axes['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = other_tiller_axes[other_tiller_axes['n0'] >= HS_break].index
        other_tiller_axes['t0'][n0_smaller_than_HS_break_indexes] = other_tiller_axes['TT_col_0'][n0_smaller_than_HS_break_indexes] + other_tiller_axes['n0'][n0_smaller_than_HS_break_indexes] / other_tiller_axes['a_cohort'][n0_smaller_than_HS_break_indexes]
        other_tiller_axes['t0'][n0_greater_than_HS_break_indexes] = (other_tiller_axes['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + other_tiller_axes['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # calculation of c
    other_tiller_axes['c'] = most_frequent_MS['c'][0] * other_tiller_axes['n2'] / most_frequent_MS['n2'][0]
    # calculation of a
    other_tiller_axes['a'] = most_frequent_MS['a'][0] * other_tiller_axes['n2'] / most_frequent_MS['n2'][0]
    # calculation of RMSE_gl
    other_tiller_axes['RMSE_gl'] = most_frequent_MS['RMSE_gl'][0]  
    return other_tiller_axes
    

def calculate_decimal_elongated_internode_number(dimT_tmp, dynT_tmp):
    '''
    Calculate the number of elongated internodes.
    
    :Parameters:
         - `dimT_tmp` (:class:`pandas.DataFrame`) - the *dimT_tmp* dataframe.
         - `dynT_tmp` (:class:`pandas.DataFrame`) - the *dynT_tmp* dataframe.

    :Returns:
        The number of elongated internodes.
    
    :Returns Type:
        :class:`float`
        
    '''
    MS_most_frequent_axis_dynT_tmp = dynT_tmp.ix[dynT_tmp.first_valid_index()]
    id_axis = MS_most_frequent_axis_dynT_tmp['id_axis']
    N_phytomer_potential = MS_most_frequent_axis_dynT_tmp['N_phytomer_potential']
    grouped = dimT_tmp.groupby(['id_axis', 'N_phytomer_potential'])
    MS_most_frequent_axis_indexes = grouped.groups[(id_axis, N_phytomer_potential)]
    # get the lengths of the internodes which belong to each phytomer of the most 
    # frequent axis of the MS
    MS_most_frequent_axis_L_internode = dimT_tmp['L_internode'][MS_most_frequent_axis_indexes]
    # keep only the non-zero lengths
    MS_most_frequent_axis_L_internode = MS_most_frequent_axis_L_internode[MS_most_frequent_axis_L_internode != 0.0]
    # get the indexes of the phytomers of the MS most frequent axis, for which L_internode is non-null
    MS_most_frequent_axis_index_phytomer = dimT_tmp['index_phytomer'][MS_most_frequent_axis_L_internode.index]
    # Fit a polynomial of degree 2 to points (MS_most_frequent_axis_L_internode, 
    # MS_most_frequent_axis_index_phytomer), and get the coefficient of degree 0. 
    return np.polyfit(MS_most_frequent_axis_L_internode, MS_most_frequent_axis_index_phytomer, 2)[2]
    
        