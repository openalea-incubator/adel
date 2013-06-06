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


def create_dynT_tmp(id_phen):
    '''
    Create the *dynT_tmp* dataframe.

    :Parameters:
    
        - `id_phen` (:class:`list` of :class:`float`) - the *id_phen* column of the :ref:`axeT <axeT>` 
          dataframe.

    :Returns:
        The *dynT_tmp* dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    ''' 
    # create id_axis from id_phen by deleting any duplicated value.
    id_axis = list(set(id_phen))
    
    # for each element in id_axis:
    # 1. decode the N_cohort
    N_cohort = [float(str(int(id_phen_))[:-2]) for id_phen_ in id_axis]
    # 2. calculate cardinality
    axis_frequencies = []
    for id_phen_ in id_axis:
        axis_frequencies.append(id_phen.count(id_phen_))    
    # 3. decode the Nff
    Nff = [float(str(int(id_phen_))[-2:]) for id_phen_ in id_axis]
    # 4. initialize the values of the other columns to NaN
    a_cohort = [np.nan for i in range(len(id_axis))]
    TT_col_0 = [np.nan for i in range(len(id_axis))]
    TT_col_break = [np.nan for i in range(len(id_axis))]
    TT_col_nff = [np.nan for i in range(len(id_axis))]
    dTT_MS_cohort = [np.nan for i in range(len(id_axis))]
    n0 = [np.nan for i in range(len(id_axis))]
    n1 = [np.nan for i in range(len(id_axis))]
    n2 = [np.nan for i in range(len(id_axis))]
    t0 = [np.nan for i in range(len(id_axis))]
    t1 = [np.nan for i in range(len(id_axis))]
    hs_t1 = [np.nan for i in range(len(id_axis))]
    a = [np.nan for i in range(len(id_axis))]
    c = [np.nan for i in range(len(id_axis))]
    RMSE_gl = [np.nan for i in range(len(id_axis))]
    
    # create the table with the appropriate columns.
    unsorted_dynT_tmp_array = np.array([N_cohort, 
                                        id_axis, 
                                        axis_frequencies, 
                                        Nff, 
                                        a_cohort, 
                                        TT_col_0, 
                                        TT_col_break, 
                                        TT_col_nff, 
                                        dTT_MS_cohort, 
                                        n0, 
                                        n1, 
                                        n2, 
                                        t0, 
                                        t1, 
                                        hs_t1, 
                                        a, 
                                        c, 
                                        RMSE_gl]).transpose()
    unsorted_dynT_tmp_dataframe = \
        pandas.DataFrame(unsorted_dynT_tmp_array, 
                          columns=['N_cohort', 
                                   'id_axis', 
                                   'cardinality', 
                                   'Nff', 
                                   'a_cohort', 
                                   'TT_col_0', 
                                   'TT_col_break', 
                                   'TT_col_nff', 
                                   'dTT_MS_cohort', 
                                   'n0', 
                                   'n1', 
                                   'n2', 
                                   't0', 
                                   't1', 
                                   'hs_t1', 
                                   'a', 
                                   'c', 
                                   'RMSE_gl'])
    # sort the table according to N_cohort (ascending order), then cardinality 
    # (descending order).
    dynT_tmp_dataframe = pandas.DataFrame(
        columns=unsorted_dynT_tmp_dataframe.columns)
    for name, group in unsorted_dynT_tmp_dataframe.groupby('N_cohort'):
        sorted_group = group.sort_index(by='cardinality', ascending=False)
        dynT_tmp_dataframe = \
            dynT_tmp_dataframe.append(sorted_group)
    # reconstruct the index  
    dynT_tmp_dataframe.index = \
        range(dynT_tmp_dataframe.index.size)
    return dynT_tmp_dataframe


def create_dynT(
       dynT_user, 
       dimT_user, 
       GL_number, 
       decimal_elongated_internode_number,
       leaf_number_delay_MS_cohort=params.LEAF_NUMBER_DELAY_MS_COHORT):
    '''
    Create the :ref:`dynT <dynT>` dataframe.
    
    :Parameters:
         - `dynT_user` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe set by the user.
         - `dimT_user` (:class:`pandas.DataFrame`) - the :ref:`dimT <dimT>` dataframe set by the user.
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
    
        * in *dynT_user*, the column *N_cohort*, *id_axis*, *cardinality* 
          and *Nff* must be completely filled, i.e. they must not contain 
          any NA value.
    
    ''' 
    # in 'dynT_user', check that the columns 'N_cohort',
    # 'id_axis', 'cardinality' and 'Nff' are non-NA.
    tools.checkValidity(dynT_user['N_cohort'].count() == dynT_user['N_cohort'].size)
    tools.checkValidity(dynT_user['id_axis'].count() == dynT_user['id_axis'].size)
    tools.checkValidity(dynT_user['cardinality'].count() == dynT_user['cardinality'].size)
    tools.checkValidity(dynT_user['Nff'].count() == dynT_user['Nff'].size)
    # get all main stem rows
    MS = dynT_user[dynT_user['N_cohort'] == 1.0]
    
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
    tiller_axes_df = dynT_user[dynT_user['N_cohort'] != 1.0]
    # extract the rows corresponding to the most frequent tiller axes
    grouped = tiller_axes_df.groupby('N_cohort')
    most_frequent_tiller_axes = []
    for N_cohort, group_indexes in grouped.groups.iteritems():
        most_frequent_tiller_axes.append(tiller_axes_df.ix[group_indexes[0:1]])
    # concatenate these rows in one dataframe ; 'most_frequent_tiller_axes' is  
    # now a pandas.DataFrame (and is not a list anymore)
    most_frequent_tiller_axes = pandas.concat(most_frequent_tiller_axes)
    # in this new dataframe, fill the columns referring to the dynamic of Haun Stage
    most_frequent_tiller_axes = _gen_most_frequent_tiller_axes_HS_dynamic(most_frequent_MS, most_frequent_tiller_axes)
    # and fill the columns referring to the dynamic of the green leaves
    most_frequent_tiller_axes = _gen_most_frequent_tiller_axes_GL_dynamic(most_frequent_MS, most_frequent_tiller_axes)
    # extract the rows corresponding to all tiller axes except the most frequent ones
    other_tiller_axes = tiller_axes_df.drop(most_frequent_tiller_axes.index)
    other_tiller_axes = _gen_other_tiller_axes_HS_dynamic(most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes)
    # and fill the columns referring to the dynamic of the green leaves
    other_tiller_axes = _gen_other_tiller_axes_GL_dynamic(most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes)
    
    return pandas.concat([most_frequent_MS, other_MS, most_frequent_tiller_axes, other_tiller_axes]).sort()
    

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
        a2_0 = (most_frequent_MS['Nff'][0] - HS_break_0) / (most_frequent_MS['TT_col_nff'][0] - most_frequent_MS['TT_col_break'][0])
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
        a2 = (most_frequent_MS['Nff'] - most_frequent_MS['HS_break']) / (most_frequent_MS['TT_col_nff'] - most_frequent_MS['TT_col_break'])
        n0_smaller_than_HS_break_indexes = most_frequent_MS[most_frequent_MS['n0'] < HS_break].index
        n0_greater_than_HS_break_indexes = most_frequent_MS[most_frequent_MS['n0'] >= HS_break].index
        most_frequent_MS['t0'][n0_smaller_than_HS_break_indexes] = most_frequent_MS['TT_col_0'][n0_smaller_than_HS_break_indexes] + most_frequent_MS['n0'][n0_smaller_than_HS_break_indexes] / most_frequent_MS['a_cohort'][n0_smaller_than_HS_break_indexes]
        most_frequent_MS['t0'][n0_greater_than_HS_break_indexes] = (most_frequent_MS['n0'][n0_greater_than_HS_break_indexes] - HS_break[n0_greater_than_HS_break_indexes]) / a2[n0_greater_than_HS_break_indexes] + most_frequent_MS['TT_col_break'][n0_greater_than_HS_break_indexes]    
    # calculation of c 
    most_frequent_MS['c'] = -((most_frequent_MS['Nff'] - decimal_elongated_internode_number) - (most_frequent_MS['n2'] - most_frequent_MS['n1'])) / (most_frequent_MS['TT_col_nff'] - most_frequent_MS['t1'])
    # calculation of a
    TT_col_nff_0 = most_frequent_MS['TT_col_nff'][0]
    n2_0 = most_frequent_MS['n2'][0]
    TT = np.array([TT_col_nff_0] + GL_number.keys()) - TT_col_nff_0
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
        other_MS['dTT_MS_cohort'] = most_frequent_MS['dTT_MS_cohort'][0] + (other_MS['id_axis'] - most_frequent_MS['id_axis'][0]) / (4 * most_frequent_MS['a_cohort'][0])
        # calculation of TT_col_nff
        other_MS['TT_col_nff'] = other_MS['dTT_MS_cohort'] + most_frequent_MS['TT_col_nff'][0]
        # calculation of a_cohort
        if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
            other_MS['a_cohort'] = other_MS['Nff'] / (other_MS['TT_col_nff'] - other_MS['TT_col_0'])
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
        other_MS['n1'] = most_frequent_MS['n1'][0] * other_MS['Nff'] / most_frequent_MS['Nff'][0]
    # calculation of t1
    other_MS['t1'] = most_frequent_MS['t1'][0] + other_MS['dTT_MS_cohort']
    # calculation of hs_t1
    other_MS['hs_t1'] = other_MS['a_cohort'] * (other_MS['t1'] - other_MS['TT_col_0'])
    # calculation of n0 
    if other_MS['n0'].count() != other_MS['n0'].size:
        other_MS['n0'] = most_frequent_MS['n0'][0] * other_MS['Nff'] / most_frequent_MS['Nff'][0]
    # calculation of n2
    if other_MS['n2'].count() != other_MS['n2'].size:
        other_MS['n2'] = most_frequent_MS['n2'][0] * other_MS['Nff'] / most_frequent_MS['Nff'][0]
    # calculation of t0
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        other_MS['t0'] = other_MS['TT_col_0'] + other_MS['n0'] / other_MS['a_cohort']
    else: # bilinear mode
        HS_break = other_MS['a_cohort'] * (other_MS['TT_col_break'] - other_MS['TT_col_0'])
        a2 = (other_MS['Nff'] - other_MS['HS_break']) / (other_MS['TT_col_nff'] - other_MS['TT_col_break'])
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
    cohorts = most_frequent_tiller_axes['N_cohort'].ix[nan_most_frequent_tiller_axis_indexes].astype(int).values
    leaf_number_delay_MS_cohorts = np.array([params.LEAF_NUMBER_DELAY_MS_COHORT[cohort] for cohort in cohorts])
    most_frequent_tiller_axes['TT_col_0'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_MS['TT_col_0'][0] + (leaf_number_delay_MS_cohorts / most_frequent_MS['a_cohort'][0])
    # dTT_MS_cohort is set by the user. Thus there is nothing to do.
    # calculation of TT_col_nff
    most_frequent_tiller_axes['TT_col_nff'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_tiller_axes['dTT_MS_cohort'].ix[nan_most_frequent_tiller_axis_indexes] + most_frequent_MS['TT_col_nff'][0]
    # calculation of a_cohort
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_tiller_axes['a_cohort'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_tiller_axes['Nff'].ix[nan_most_frequent_tiller_axis_indexes] / (most_frequent_tiller_axes['TT_col_nff'].ix[nan_most_frequent_tiller_axis_indexes] - most_frequent_tiller_axes['TT_col_0'].ix[nan_most_frequent_tiller_axis_indexes])
    else: # bilinear mode
        most_frequent_tiller_axes['a_cohort'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_MS['a_cohort'][0]
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
    most_frequent_tiller_axes['n1'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_MS['n1'][0]
    # calculation of t1
    most_frequent_tiller_axes['t1'] = most_frequent_MS['t1'][0] + most_frequent_tiller_axes['dTT_MS_cohort']
    # calculation of hs_t1
    most_frequent_tiller_axes['hs_t1'] = most_frequent_tiller_axes['a_cohort'] * (most_frequent_tiller_axes['t1'] - most_frequent_tiller_axes['TT_col_0'])
    # calculation of n0
    n1_hs_t1 = most_frequent_tiller_axes[['n1', 'hs_t1']].ix[nan_most_frequent_tiller_axis_indexes]
    if n1_hs_t1.index.size != 0:
        most_frequent_tiller_axes['n0'].ix[nan_most_frequent_tiller_axis_indexes] = n1_hs_t1.apply(np.min, 1)                      
    # calculation of n2
    most_frequent_tiller_axes['n2'].ix[nan_most_frequent_tiller_axis_indexes] = most_frequent_MS['n2'][0] * params.N2_MS_DIV_N2_COHORT
    # calculation of t0
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        most_frequent_tiller_axes['t0'] = most_frequent_tiller_axes['TT_col_0'] + most_frequent_tiller_axes['n0'] / most_frequent_tiller_axes['a_cohort']
    else: # bilinear mode
        HS_break = most_frequent_tiller_axes['a_cohort'] * (most_frequent_tiller_axes['TT_col_break'] - most_frequent_tiller_axes['TT_col_0'])
        a2 = (most_frequent_tiller_axes['Nff'] - most_frequent_tiller_axes['HS_break']) / (most_frequent_tiller_axes['TT_col_nff'] - most_frequent_tiller_axes['TT_col_break'])
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
    for name, group in other_tiller_axes.ix[nan_other_tiller_axis_indexes].groupby('N_cohort'):
        most_frequent_tiller_axis_idx = most_frequent_tiller_axes.ix[most_frequent_tiller_axes['N_cohort'] == name].first_valid_index()
        # calculation of TT_col_0
        other_tiller_axes['TT_col_0'].ix[group.index] = most_frequent_tiller_axes['TT_col_0'][most_frequent_tiller_axis_idx]
        # calculation of dTT_MS_cohort
        other_tiller_axes['dTT_MS_cohort'].ix[group.index] = \
            most_frequent_tiller_axes['dTT_MS_cohort'][most_frequent_tiller_axis_idx] \
            + (group['id_axis'] - most_frequent_tiller_axes['id_axis'][most_frequent_tiller_axis_idx]) \
            / (4 * most_frequent_MS['a_cohort'][0])
        if most_frequent_MS['TT_col_break'][0] != 0.0: 
            # calculation of a_cohort in bilinear mode
            other_tiller_axes['a_cohort'].ix[group.index] = most_frequent_tiller_axes['a_cohort'][most_frequent_tiller_axis_idx]
    # calculation of TT_col_nff
    other_tiller_axes['TT_col_nff'].ix[nan_other_tiller_axis_indexes] = other_tiller_axes['dTT_MS_cohort'].ix[nan_other_tiller_axis_indexes] + most_frequent_MS['TT_col_nff'][0]
    # calculation of a_cohort in linear mode
    if most_frequent_MS['TT_col_break'][0] == 0.0:
        other_tiller_axes['a_cohort'].ix[nan_other_tiller_axis_indexes] = other_tiller_axes['Nff'].ix[nan_other_tiller_axis_indexes] / (other_tiller_axes['TT_col_nff'].ix[nan_other_tiller_axis_indexes] - other_tiller_axes['TT_col_0'].ix[nan_other_tiller_axis_indexes])
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
    for name, group in other_tiller_axes.ix[nan_other_tiller_axis_indexes].groupby('N_cohort'):
        most_frequent_tiller_axis_idx = most_frequent_tiller_axes.ix[most_frequent_tiller_axes['N_cohort'] == name].first_valid_index()
        # calculation ofn1
        other_tiller_axes['n1'].ix[group.index] = most_frequent_tiller_axes['n1'][most_frequent_tiller_axis_idx]
        # calculation of n2
        other_tiller_axes['n2'].ix[group.index] = most_frequent_tiller_axes['n2'][most_frequent_tiller_axis_idx]
    # calculation of t1
    other_tiller_axes['t1'] = most_frequent_MS['t1'][0] + other_tiller_axes['dTT_MS_cohort']
    # calculation of hs_t1
    other_tiller_axes['hs_t1'] = other_tiller_axes['a_cohort'] * (other_tiller_axes['t1'] - other_tiller_axes['TT_col_0'])
    # calculation of n0
    n1_hs_t1 = other_tiller_axes[['n1', 'hs_t1']].ix[nan_other_tiller_axis_indexes]
    if n1_hs_t1.index.size != 0:
        other_tiller_axes['n0'].ix[nan_other_tiller_axis_indexes] = n1_hs_t1.apply(np.min, 1)                      
    # calculation of t0
    if most_frequent_MS['TT_col_break'][0] == 0.0: # linear mode
        other_tiller_axes['t0'] = other_tiller_axes['TT_col_0'] + other_tiller_axes['n0'] / other_tiller_axes['a_cohort']
    else: # bilinear mode
        HS_break = other_tiller_axes['a_cohort'] * (other_tiller_axes['TT_col_break'] - other_tiller_axes['TT_col_0'])
        a2 = (other_tiller_axes['Nff'] - other_tiller_axes['HS_break']) / (other_tiller_axes['TT_col_nff'] - other_tiller_axes['TT_col_break'])
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
    

def calculate_decimal_elongated_internode_number(dimT_user):
    '''
    Calculate the number of elongated internodes.
    
    :Parameters:
         - `dimT_user` (:class:`pandas.DataFrame`) - the :ref:`dimT <dimT>` 
           dataframe set by the user.

    :Returns:
        The number of elongated internodes.
    
    :Returns Type:
        :class:`float`
        
    .. warning:: *dimT_user* must be a :class:`pandas.DataFrame`.
    
    '''
    # calculate the number of rows which describe the first axis 
    first_axis_rows_number = int(str(int(dimT_user['id_dim'][0]))[-2:])
    # get the length of the internodes which belong to each phytomer of the first axis
    first_axis_L_internode = dimT_user['L_internode'][:first_axis_rows_number]
    # keep only the non-zero lengths
    first_axis_L_internode = first_axis_L_internode[first_axis_L_internode != 0.0]
    # get the index of each phytomer which belongs to the first axis  
    first_axis_index_phytomer = dimT_user['index_phytomer'][first_axis_L_internode.index]
    # Fit a polynomial of degree 2 to points (first_axis_L_internode, 
    # first_axis_index_phytomer), and get the coefficient of degree 0. 
    return np.polyfit(first_axis_L_internode, first_axis_index_phytomer, 2)[2]
    
    