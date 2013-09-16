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
Provides functions to construct the *dimT_tmp*, the :ref:`dimT_abs <dimT_abs>` and the :ref:`dimT <dimT>` dataframes.

The :ref:`dimT_abs <dimT_abs>` and the :ref:`dimT <dimT>` dataframes are described in the User Guide 
(see :ref:`adel_user`).

Authors: Mariem Abichou, Camille Chambon, Bruno Andrieu
'''

import pandas
import numpy as np

from adel.plantgen import tools, params


def create_dimT_tmp(axeT_tmp):
    '''
    Create the *dimT_tmp* dataframe.
    Compute the following columns: *id_axis*, *N_phytomer_potential*, *index_phytomer*.
    
    :Parameters:
    
        - `axeT_tmp` (:class:`pandas.DataFrame`) - the *axeT_tmp* dataframe.
        
    :Returns: 
        The *dimT_tmp* dataframe.
    
    :Returns Type: 
        :class:`pandas.DataFrame`
    
    '''
    id_axis_list = []
    N_phytomer_potential_list = []
    index_phytomer_list = []
    for (id_axis, N_phytomer_potential), axeT_tmp_group in axeT_tmp.groupby(['id_axis', 'N_phytomer_potential']):
        id_axis_list.extend(np.repeat(id_axis, N_phytomer_potential))
        N_phytomer_potential_list.extend(np.repeat(N_phytomer_potential, N_phytomer_potential))
        index_phytomer_list.extend(range(1, N_phytomer_potential + 1))
    
    dimT_tmp = pandas.DataFrame(index=range(len(id_axis_list)),
                                columns=['id_axis', 'N_phytomer_potential', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'],
                                dtype=float)
    dimT_tmp['id_axis'] = id_axis_list
    dimT_tmp['N_phytomer_potential'] = N_phytomer_potential_list
    dimT_tmp['index_phytomer'] = index_phytomer_list
    
    return dimT_tmp
    

def create_dimT_abs(axeT_, dimT_tmp, phenT_abs, dynT_):
    '''
    Create the :ref:`dimT_abs <dimT_abs>` dataframe filling the *dimT_tmp* dataframe.

    :Parameters:
        
        - `axeT_` (:class:`pandas.DataFrame`) - the :ref:`axeT <axeT>` dataframe.
        - `dimT_tmp` (:class:`pandas.DataFrame`) - the *dimT_tmp* dataframe.
        - `phenT_abs` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
        - `dynT_` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe.
    
    :Returns: 
        The :ref:`dimT_abs <dimT_abs>` dataframe.
    
    :Returns Type: 
        :class:`pandas.DataFrame`
        
    .. warning:: 
        * in *dimT_tmp*, the columns *id_axis*, *N_phytomer_potential* and *index_phytomer* 
          must be completely filled, i.e. they must not contain any NA value.
        * in *dimT_tmp*, the rows which describe the most frequent axis of the MS 
          must be completely filled, i.e. there must not contain any NA value.
        
    '''
    if dimT_tmp['id_axis'].count() != dimT_tmp['id_axis'].size:
        raise tools.InputError("dimT_tmp['id_axis'] contains NA values")
    if dimT_tmp['N_phytomer_potential'].count() != dimT_tmp['N_phytomer_potential'].size:
        raise tools.InputError("dimT_tmp['N_phytomer_potential'] contains NA values")
    if dimT_tmp['index_phytomer'].count() != dimT_tmp['index_phytomer'].size:
        raise tools.InputError("dimT_tmp['index_phytomer'] contains NA values")
    
    dimT_tmp_grouped = dimT_tmp.groupby(['id_axis', 'N_phytomer_potential'])
    dimT_tmp_group = dimT_tmp_grouped.get_group((dynT_['id_axis'][0], dynT_['N_phytomer_potential'][0]))
    dimT_tmp_group_without_na = dimT_tmp_group.dropna()
    if len(dimT_tmp_group_without_na) != len(dimT_tmp_group):
        raise tools.InputError("dimT_tmp does not contain the dimensions of the most frequent MS")
    
    dimT_abs = _init_dimT_abs(axeT_, 
                              dimT_tmp, 
                              phenT_abs, 
                              dynT_)
    
    MS_dynT = dynT_[dynT_['id_axis'] == 'MS']
    idxmax = MS_dynT['cardinality'].idxmax()
    MS_id_cohort = MS_dynT['id_cohort'][idxmax]
    MS_N_phytomer_potential = MS_dynT['N_phytomer_potential'][idxmax]
    axeT_grouped = axeT_.groupby(['id_axis', 'id_cohort', 'N_phytomer_potential'])
    axeT_group = axeT_grouped.get_group(('MS', MS_id_cohort, MS_N_phytomer_potential))
    MS_id_dim = axeT_group['id_dim'][axeT_group.first_valid_index()]

    L_blade_is_null = dimT_abs['L_blade'].isnull()
    row_indexes_to_fit = L_blade_is_null[L_blade_is_null == True].index
    
    _gen_lengths(MS_id_dim, row_indexes_to_fit, dimT_abs)
    
    _gen_widths(MS_id_dim, row_indexes_to_fit, dimT_abs)
    
    dimT_abs.sort(['is_ear', 'id_dim'], inplace=True)
    
    # reinitialize the index
    dimT_abs.index = range(dimT_abs.index.size)
            
    return dimT_abs.drop(['TT_em_phytomer', 'is_ear'], axis=1)


def _init_dimT_abs(axeT_, dimT_tmp, phenT_abs, dynT_):
    '''Initialize dimT_abs.'''
    # create a new dataframe from phenT_abs, removing the lines for which index_phytomer==0.0, 
    # and keeping only the columns 'id_phen', 'index_phytomer' and 'TT_em_phytomer'.
    TT_em_phytomer_cleaned = phenT_abs.drop(phenT_abs.groupby('index_phytomer').groups[0.0])['TT_em_phytomer']
    id_phen_cleaned = phenT_abs['id_phen'][TT_em_phytomer_cleaned.index].astype(int)
    index_phytomer_cleaned = phenT_abs['index_phytomer'][TT_em_phytomer_cleaned.index].astype(int)
    phenT_abs_cleaned = pandas.DataFrame(np.array([id_phen_cleaned.values, index_phytomer_cleaned.values, TT_em_phytomer_cleaned.values]).transpose(),
                                         columns=['id_phen', 'index_phytomer', 'TT_em_phytomer'])
    
    dimT_abs = pandas.DataFrame(columns=['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode', 'TT_em_phytomer', 'is_ear'])
    
    phenT_abs_cleaned_grouped = phenT_abs_cleaned.groupby('id_phen')
    dimT_tmp_grouped = dimT_tmp.groupby(['id_axis', 'N_phytomer_potential'])
    for id_dim, axeT_group in axeT_.groupby('id_dim'):
        axeT_keys = axeT_group.groupby(['id_axis', 'id_cohort', 'N_phytomer_potential']).groups.keys()
        dynT_group = dynT_.select(lambda idx: (dynT_['id_axis'][idx], dynT_['id_cohort'][idx], dynT_['N_phytomer_potential'][idx]) in axeT_keys)
        idxmax = dynT_group[dynT_group['id_axis'] == dynT_group['id_axis'].max()].first_valid_index()
        id_cohort = dynT_group['id_cohort'][idxmax]
        N_phytomer_potential = dynT_group['N_phytomer_potential'][idxmax]
        id_phen = int(''.join([str(id_cohort), str(N_phytomer_potential).zfill(2)]))
            
        phenT_abs_cleaned_group = phenT_abs_cleaned_grouped.get_group(id_phen)
        
        dimT_abs_group_idx = np.arange(axeT_group['N_phytomer'][axeT_group.first_valid_index()])
        
        dimT_abs_group = pandas.DataFrame(index=dimT_abs_group_idx, 
                                          columns=dimT_abs.columns,
                                          dtype=float)
        dimT_abs_group['id_dim'] = id_dim
        dimT_abs_group['index_phytomer'] = dimT_abs_group_idx + 1
        
        is_ear = int(str(int(id_dim))[-1])
        
        if is_ear == 1:
            id_axis = dynT_group['id_axis'][idxmax]
            dimT_tmp_group = dimT_tmp_grouped.get_group((id_axis, N_phytomer_potential))
            organ_dim_list = ['L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
            for organ_dim in organ_dim_list:
                dim_idx_to_get = dimT_tmp_group.index[dimT_abs_group.index]
                dimT_abs_group[organ_dim] = dimT_tmp_group[organ_dim][dim_idx_to_get].values.astype(float)
        
        phen_idx_to_get = phenT_abs_cleaned_group.index[dimT_abs_group.index]
        dimT_abs_group['TT_em_phytomer'] = phenT_abs_cleaned_group['TT_em_phytomer'][phen_idx_to_get].values.astype(float)
        
        dimT_abs_group['is_ear'] = is_ear
        
        dimT_abs = pandas.concat([dimT_abs, dimT_abs_group], ignore_index=True)
    
    # force the type of id_dim and is_ear
    dimT_abs[['id_dim', 'index_phytomer', 'is_ear']] = dimT_abs[['id_dim', 'index_phytomer', 'is_ear']].astype(int)

    return dimT_abs
    

def _gen_lengths(MS_id_dim, row_indexes_to_fit, dimT_abs):
    '''Fit the lengths in-place.'''
    TT_em_phytomer_series = dimT_abs['TT_em_phytomer']
    MS_rows_indexes = dimT_abs[dimT_abs['id_dim'] == MS_id_dim].index
    MS_last_TT_em_phytomer = TT_em_phytomer_series[MS_rows_indexes[-1]]
    indexes_to_ceil = TT_em_phytomer_series[TT_em_phytomer_series > MS_last_TT_em_phytomer].index
    
    for length in ['L_blade', 'L_sheath', 'L_internode']:
        current_length_series = dimT_abs[length]
        MS_lengths_series = current_length_series[MS_rows_indexes]
        MS_non_null_lengths_rows_indexes = MS_lengths_series[MS_lengths_series != 0.0].index
        polynomial_coefficients_array = np.polyfit(dimT_abs['TT_em_phytomer'][MS_non_null_lengths_rows_indexes].values, 
                                                   current_length_series[MS_non_null_lengths_rows_indexes].values, 6)
        MS_last_length = MS_lengths_series[MS_non_null_lengths_rows_indexes[-1]]
        for id_dim, dimT_abs_group in dimT_abs.ix[row_indexes_to_fit].groupby(by='id_dim'):
            TT_em_phytomer_group = dimT_abs_group['TT_em_phytomer']
            current_length_series[dimT_abs_group.index] = np.polyval(polynomial_coefficients_array, 
                                                                     TT_em_phytomer_group)
            # ceiling
            current_length_series[indexes_to_ceil.intersection(dimT_abs_group.index)] = MS_last_length
            if dimT_abs_group['is_ear'][dimT_abs_group.first_valid_index()] == 0:
                # regression
                for i in range(len(params.REGRESSION_OF_DIMENSIONS[length])):
                    current_phytomer_index = i-3
                    if current_phytomer_index in dimT_abs_group.index:
                        current_length_series[dimT_abs_group.index[current_phytomer_index]] *= (1.0 - params.REGRESSION_OF_DIMENSIONS[length][i])

        # thresholding
        if length == 'L_internode':
            MS_first_non_null_TT_em_phytomer = TT_em_phytomer_series[MS_non_null_lengths_rows_indexes[0]]
            indexes_to_threshold = TT_em_phytomer_series[TT_em_phytomer_series <= MS_first_non_null_TT_em_phytomer].index
            indexes_to_threshold = indexes_to_threshold.intersection(row_indexes_to_fit)
            current_length_series[indexes_to_threshold] = 0.0


def _gen_widths(MS_id_dim, row_indexes_to_fit, dimT_abs):
    '''Fit the widths in-place.'''
    MS_rows_indexes = dimT_abs[dimT_abs['id_dim'] == MS_id_dim].index
    TT_em_phytomer_series = dimT_abs['TT_em_phytomer']
    for width in ['W_blade', 'W_sheath', 'W_internode']:
        current_width_series = dimT_abs[width]
        MS_width_series = current_width_series[MS_rows_indexes]
        MS_non_null_widths_rows_indexes = MS_width_series[MS_width_series != 0.0].index
        MS_first_non_null_width = current_width_series[MS_non_null_widths_rows_indexes[0]]
        MS_last_non_null_width = current_width_series[MS_non_null_widths_rows_indexes[-1]]
        MS_first_non_null_TT_em_phytomer = TT_em_phytomer_series[MS_non_null_widths_rows_indexes[0]]
        MS_last_non_null_TT_em_phytomer = TT_em_phytomer_series[MS_non_null_widths_rows_indexes[-1]]
        
        for id_dim, dimT_abs_group in dimT_abs.ix[row_indexes_to_fit].groupby(by='id_dim'):
            TT_em_phytomer_group = dimT_abs_group['TT_em_phytomer']
            if width == 'W_internode':
                # get TT_em_phytomer of the main stem first phytomer which has a 
                # TT_em_phytomer greater than the first main stem phytomer with 
                # a non null width
                valid_TT_em_phytomers = TT_em_phytomer_group[TT_em_phytomer_group >= MS_first_non_null_TT_em_phytomer]
                if len(valid_TT_em_phytomers) == 0:
                    continue # The widths of these phytomers are thresholded to 0 later.
                x1 = valid_TT_em_phytomers[valid_TT_em_phytomers.index[0]]
                # get TT_em_phytomer of the main stem last phytomer which has a 
                # TT_em_phytomer lesser than the last main stem phytomer
                valid_TT_em_phytomers = TT_em_phytomer_group[TT_em_phytomer_group <= MS_last_non_null_TT_em_phytomer]
                x2 = valid_TT_em_phytomers[valid_TT_em_phytomers.index[-1]]
            else:
                x1 = TT_em_phytomer_group[TT_em_phytomer_group.first_valid_index()]
                x2 = TT_em_phytomer_group[TT_em_phytomer_group.last_valid_index()]
                
            y1 = MS_first_non_null_width
            y2 = MS_last_non_null_width
            polynomial_coefficient_array = np.polyfit(np.array([x1, x2]), np.array([y1, y2]), 1)
            current_width_series[dimT_abs_group.index] = np.polyval(polynomial_coefficient_array, TT_em_phytomer_group)
            if width == 'W_internode':
                # ceiling of the width of the phytomers which have a TT_em_phytomer 
                # greater than the TT_em_phytomer of the last main stem 
                indexes_to_ceil = current_width_series[current_width_series > MS_last_non_null_TT_em_phytomer].index
                current_width_series[indexes_to_ceil.intersection(dimT_abs_group.index)] = MS_last_non_null_width
            
            if dimT_abs_group['is_ear'][dimT_abs_group.first_valid_index()] == 0:
                # regression
                for i in range(len(params.REGRESSION_OF_DIMENSIONS[width])):
                    current_phytomer_index = i-3
                    if current_phytomer_index in dimT_abs_group.index:
                        current_width_series[dimT_abs_group.index[current_phytomer_index]] *= (1.0 - params.REGRESSION_OF_DIMENSIONS[width][i])
        
        if width == 'W_internode':
            # thresholding of the width of the phytomers which have a TT_em_phytomer 
            # lesser than the TT_em_phytomer of the first main stem which has a non null width
            indexes_to_threshold = TT_em_phytomer_series[TT_em_phytomer_series <= MS_first_non_null_TT_em_phytomer].index
            indexes_to_threshold = indexes_to_threshold.intersection(row_indexes_to_fit)
            current_width_series[indexes_to_threshold] = 0.0
            
        current_width_series = current_width_series.clip_lower(0.0)
    

def create_dimT(dimT_abs):
    '''
    Create the :ref:`dimT <dimT>` dataframe. 
    
    :Parameters:
    
        - `dimT_abs` (:class:`pandas.DataFrame`) - the :ref:`dimT_abs <dimT_abs>` dataframe.
        
    :Returns: 
        the :ref:`dimT <dimT>` dataframe.
        
    :Returns Type: 
        pandas.Dataframe
        
    .. warning:: 
        
        * *dimT_abs* must be completely filled, i.e. must not contain any 
          NA value.
                 
    '''
    if not (dimT_abs.count().max() == dimT_abs.count().min() == dimT_abs.index.size):
        raise tools.InputError("dimT_abs contains NA values")

    dimT_ = pandas.DataFrame(index=dimT_abs.index, columns=['id_dim', 'index_rel_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'], dtype=float)
    dimT_[['id_dim', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']] = dimT_abs[['id_dim', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']]
    tmp_series = pandas.Series(dimT_.index)
    for name, group in dimT_abs.groupby('id_dim'):
        def normalize_index_phytomer(i):
            return group['index_phytomer'][i] / float(group['index_phytomer'][group.index[-1]])
        dimT_['index_rel_phytomer'][group.index] = tmp_series[group.index].map(normalize_index_phytomer) 
        
    return dimT_
        
        
    