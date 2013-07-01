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


def create_dimT_tmp(dynT_tmp):
    '''
    Create the *dimT_tmp* dataframe.
    Compute the following columns: *id_dim*, *index_phytomer*.
    
    :Parameters:
    
        - `dynT_tmp` (:class:`pandas.DataFrame`) - the *dynT_tmp* dataframe.
        
    :Returns: 
        The *dimT_tmp* dataframe.
    
    :Returns Type: 
        :class:`pandas.DataFrame`
    
    '''
    id_dim_tmp_list = _gen_id_dim_tmp_list(dynT_tmp)
    index_phytomer_tmp_list = _gen_index_phytomer_tmp_list(id_dim_tmp_list)
    L_blade_list = [np.nan for i in range(len(id_dim_tmp_list))]
    W_blade_list = [np.nan for i in range(len(id_dim_tmp_list))]
    L_sheath_list = [np.nan for i in range(len(id_dim_tmp_list))]
    W_sheath_list = [np.nan for i in range(len(id_dim_tmp_list))]
    L_internode_list = [np.nan for i in range(len(id_dim_tmp_list))]
    W_internode_list = [np.nan for i in range(len(id_dim_tmp_list))]
    dimT_tmp_array = np.array([id_dim_tmp_list, index_phytomer_tmp_list, L_blade_list, W_blade_list, L_sheath_list, W_sheath_list, L_internode_list, W_internode_list]).transpose()
    return pandas.DataFrame(dimT_tmp_array, columns=['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'])
    

def create_dimT_abs(axeT, dimT_user, phenT_abs, dynT):
    '''
    Create the :ref:`dimT_abs <dimT_abs>` dataframe filling the *dimT_user* dataframe.

    :Parameters:
        
        - `axeT` (:class:`pandas.DataFrame`) - the :ref:`axeT <axeT>` dataframe.
        - `dimT_user` (:class:`pandas.DataFrame`) - the :ref:`dimT <dimT>` dataframe set by the user.
        - `phenT_abs` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
        - `dynT` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe.
    
    :Returns: 
        The :ref:`dimT_abs <dimT_abs>` dataframe.
    
    :Returns Type: 
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * in *dimT_user*, the column *id_dim* and *index_phytomer* must be 
          completely filled, i.e. they must not contain any NA value.
        * in *dimT_user*, the rows which describe the first axis must be 
          completely filled, i.e. there must not contain any NA value.
        
    '''
    tools.checkValidity(dimT_user['id_dim'].count() == dimT_user['id_dim'].size)
    tools.checkValidity(dimT_user['index_phytomer'].count() == dimT_user['index_phytomer'].size)
    MS_id_axis = dynT['id_axis'][0]
    MS_row_indexes = dimT_user[dimT_user['id_dim'] == MS_id_axis].index
    for column_name in dimT_user:
        tools.checkValidity(dimT_user[column_name][MS_row_indexes].map(lambda x: x == 0.0 and 1.0 or x).fillna(0.0).all())
           
    dimT_abs_dataframe = _init_dimT_abs(axeT['id_dim'].unique(), 
                                        dimT_user, 
                                        phenT_abs, 
                                        dynT)
    
    L_blade_is_null = dimT_abs_dataframe['L_blade'].isnull()
    row_indexes_to_fit = L_blade_is_null[L_blade_is_null == True].index
    
    _gen_lengths(MS_id_axis, row_indexes_to_fit, dimT_abs_dataframe)
    
    _gen_widths(MS_id_axis, row_indexes_to_fit, dimT_abs_dataframe)
            
    return dimT_abs_dataframe.drop(['TT_em_phytomer', 'is_ear'], axis=1)


def _gen_id_dim_tmp_list(dyn_tmp_dataframe):
    '''Generate the *id_dim* column of the dimT_tmp table'''
    sorted_id_axis = dyn_tmp_dataframe['id_axis']
    id_dim_tmp_list = []
    for id_dim in sorted_id_axis:
        N_phyt = int(str(int(id_dim))[-2:])
        for i in range(N_phyt):
            id_dim_tmp_list.append(id_dim)
    return id_dim_tmp_list


def _init_dimT_abs(id_dim_from_axeT_unique, dimT_user, phenT_abs, dynT):
    '''Initialize dimT_abs.'''
    
    # create a new dateframe from phenT_abs, removing the lines for which index_phytomer==0.0, 
    # and keeping only the columns 'id_phen', 'index_phytomer' and 'TT_em_phytomer'.
    TT_em_phytomer_cleaned = phenT_abs.drop(phenT_abs.groupby('index_phytomer').groups[0.0])['TT_em_phytomer']
    id_phen_cleaned = phenT_abs['id_phen'][TT_em_phytomer_cleaned.index].astype(int)
    index_phytomer_cleaned = phenT_abs['index_phytomer'][TT_em_phytomer_cleaned.index].astype(int)
    phenT_abs_cleaned = pandas.DataFrame(np.array([id_phen_cleaned.values, index_phytomer_cleaned.values, TT_em_phytomer_cleaned.values]).transpose(),
                                         columns=['id_phen', 'index_phytomer', 'TT_em_phytomer'])
    
    # initialize dimT_abs from dimT_user
    dimT_abs = dimT_user.copy()
    dimT_abs['id_dim'] = dimT_abs['id_dim'].astype(int) * 10 + 1
    dimT_abs['TT_em_phytomer'] = phenT_abs_cleaned['TT_em_phytomer']
    dimT_abs['is_ear'] = 1
    
    # initialize dimT_abs from axeT
    id_dim_from_axeT_unique = id_dim_from_axeT_unique.astype(int)
    id_dim_from_dimT_unique = dimT_abs['id_dim'].unique()
    
    dynT_grouped = dynT.groupby('N_cohort')
    for id_dim_from_axeT_ in id_dim_from_axeT_unique:
        if id_dim_from_axeT_ not in id_dim_from_dimT_unique:
            current_id_dim_from_axeT_str = str(id_dim_from_axeT_)
            current_N_phytomer = int(current_id_dim_from_axeT_str[1:-1])
            # TODO: remove this workaround
            if current_N_phytomer == 0:
                continue
            
            id_phen_from_id_dim = int(current_id_dim_from_axeT_str[:-1])
            if id_phen_from_id_dim in phenT_abs_cleaned['id_phen'].values:
                # TT_em_phytomer of the corresponding ear.
                TT_em_phytomer_to_extract_indexes = phenT_abs_cleaned[phenT_abs_cleaned['id_phen']==id_phen_from_id_dim].index
                new_TT_em_phytomer = phenT_abs_cleaned['TT_em_phytomer'][TT_em_phytomer_to_extract_indexes].values
                index_phytomer_list = phenT_abs_cleaned['index_phytomer'][TT_em_phytomer_to_extract_indexes].tolist()
                
            else:
                # TT_em_phytomer of the most frequent axis of the current cohort.
                current_N_cohort = int(current_id_dim_from_axeT_str[:-3])
                dynT_group = dynT_grouped.get_group(current_N_cohort)
                most_frequent_axis = dynT_group['id_axis'][dynT_group.first_valid_index()]
                TT_em_phytomer_of_most_frequent_axis = phenT_abs_cleaned[phenT_abs_cleaned['id_phen']==most_frequent_axis]['TT_em_phytomer']
                number_of_TT_em_phytomer_to_extract = min(current_N_phytomer, len(TT_em_phytomer_of_most_frequent_axis))
                index_phytomer_list = range(1, number_of_TT_em_phytomer_to_extract + 1)
                TT_em_phytomer_to_extract_indexes = TT_em_phytomer_of_most_frequent_axis.index[:number_of_TT_em_phytomer_to_extract]
                new_TT_em_phytomer = TT_em_phytomer_of_most_frequent_axis[TT_em_phytomer_to_extract_indexes].values
            
            new_dimT = pandas.DataFrame(index=index_phytomer_list, 
                                        columns=dimT_abs.columns)
            new_dimT[['L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']] = np.nan
            new_dimT = new_dimT.astype(float)
            new_dimT['id_dim'] = id_dim_from_axeT_
            new_dimT['index_phytomer'] = index_phytomer_list
            new_dimT['TT_em_phytomer'] = new_TT_em_phytomer
            new_dimT['is_ear'] = int(current_id_dim_from_axeT_str[-1])
            dimT_abs = pandas.concat([dimT_abs, new_dimT], ignore_index=True)
        
    dimT_abs = dimT_abs[['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode', 'TT_em_phytomer', 'is_ear']]

    return dimT_abs
    

def _gen_index_phytomer_tmp_list(id_dim_tmp_list):
    '''Generate the *index_phytomer* column of dimT_abs.'''
    index_phytomer_tmp_list = []
    i = 0
    while i < len(id_dim_tmp_list):
        N_phyt = int(str(int(id_dim_tmp_list[i]))[-2:])
        new_index_phytomers = range(1, N_phyt + 1)
        index_phytomer_tmp_list.extend(new_index_phytomers)
        i = i + len(new_index_phytomers)
    return index_phytomer_tmp_list


def _gen_lengths(MS_id_axis, row_indexes_to_fit, dimT_abs_dataframe):
    '''Fit the lengths in-place.'''
    TT_em_phytomer_series = dimT_abs_dataframe['TT_em_phytomer']
    old_id_dim = MS_id_axis * 10 + 1
    MS_rows_indexes = dimT_abs_dataframe[dimT_abs_dataframe['id_dim'] == old_id_dim].index
    MS_last_TT_em_phytomer = TT_em_phytomer_series[MS_rows_indexes[-1]]
    indexes_to_ceil = TT_em_phytomer_series[TT_em_phytomer_series > MS_last_TT_em_phytomer].index
    
    for length in ['L_blade', 'L_sheath', 'L_internode']:
        current_length_series = dimT_abs_dataframe[length]
        MS_lengths_series = current_length_series[MS_rows_indexes]
        MS_non_null_lengths_rows_indexes = MS_lengths_series[MS_lengths_series != 0.0].index
        polynomial_coefficients_array = np.polyfit(dimT_abs_dataframe['TT_em_phytomer'][MS_non_null_lengths_rows_indexes].values, 
                                                   current_length_series[MS_non_null_lengths_rows_indexes].values, 6)
        MS_last_length = MS_lengths_series[MS_non_null_lengths_rows_indexes[-1]]
        for name, group in dimT_abs_dataframe.ix[row_indexes_to_fit].groupby(by='id_dim'):
            TT_em_phytomer_group = group['TT_em_phytomer']
            current_length_series[group.index] = np.polyval(polynomial_coefficients_array, 
                                                           TT_em_phytomer_group)
            # ceiling
            current_length_series[indexes_to_ceil.intersection(group.index)] = MS_last_length
            if group['is_ear'][group.first_valid_index()] == 0:
                # regression
                for i in range(len(params.REGRESSION_OF_DIMENSIONS[length])):
                    current_phytomer_index = i-3
                    if current_phytomer_index in group.index:
                        current_length_series[group.index[current_phytomer_index]] *= (1.0 - params.REGRESSION_OF_DIMENSIONS[length][i])

        # thresholding
        if length == 'L_internode':
            MS_first_non_null_TT_em_phytomer = TT_em_phytomer_series[MS_non_null_lengths_rows_indexes[0]]
            indexes_to_threshold = TT_em_phytomer_series[TT_em_phytomer_series <= MS_first_non_null_TT_em_phytomer].index
            indexes_to_threshold = indexes_to_threshold.intersection(row_indexes_to_fit)
            current_length_series[indexes_to_threshold] = 0.0


def _gen_widths(MS_id_axis, row_indexes_to_fit, dimT_abs_dataframe):
    '''Fit the widths in-place.'''
    old_id_dim = MS_id_axis * 10 + 1
    MS_rows_indexes = dimT_abs_dataframe[dimT_abs_dataframe['id_dim'] == old_id_dim].index
    TT_em_phytomer_series = dimT_abs_dataframe['TT_em_phytomer']
    for width in ['W_blade', 'W_sheath', 'W_internode']:
        current_width_series = dimT_abs_dataframe[width]
        MS_width_series = current_width_series[MS_rows_indexes]
        MS_non_null_widths_rows_indexes = MS_width_series[MS_width_series != 0.0].index
        MS_first_non_null_width = current_width_series[MS_non_null_widths_rows_indexes[0]]
        MS_last_non_null_width = current_width_series[MS_non_null_widths_rows_indexes[-1]]
        MS_first_non_null_TT_em_phytomer = TT_em_phytomer_series[MS_non_null_widths_rows_indexes[0]]
        MS_last_non_null_TT_em_phytomer = TT_em_phytomer_series[MS_non_null_widths_rows_indexes[-1]]
        
        for name, group in dimT_abs_dataframe.ix[row_indexes_to_fit].groupby(by='id_dim'):
            TT_em_phytomer_group = group['TT_em_phytomer']
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
            current_width_series[group.index] = np.polyval(polynomial_coefficient_array, TT_em_phytomer_group)
            if width == 'W_internode':
                # ceiling of the width of the phytomers which have a TT_em_phytomer 
                # greater than the TT_em_phytomer of the last main stem 
                indexes_to_ceil = current_width_series[current_width_series > MS_last_non_null_TT_em_phytomer].index
                current_width_series[indexes_to_ceil.intersection(group.index)] = MS_last_non_null_width
            
            if group['is_ear'][group.first_valid_index()] == 0:
                # regression
                for i in range(len(params.REGRESSION_OF_DIMENSIONS[width])):
                    current_phytomer_index = i-3
                    if current_phytomer_index in group.index:
                        current_width_series[group.index[current_phytomer_index]] *= (1.0 - params.REGRESSION_OF_DIMENSIONS[width][i])
        
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
    tools.checkValidity(dimT_abs.count().max() == dimT_abs.count().min() == dimT_abs.index.size)
    dimT_dataframe = pandas.DataFrame(index=dimT_abs.index, columns=['id_dim', 'index_rel_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'], dtype=float)
    dimT_dataframe[['id_dim', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']] = dimT_abs[['id_dim', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']]
    tmp_series = pandas.Series(dimT_dataframe.index)
    for name, group in dimT_abs.groupby('id_dim'):
        def normalize_index_phytomer(i):
            return group['index_phytomer'][i] / float(group['index_phytomer'][group.index[-1]])
        dimT_dataframe['index_rel_phytomer'][group.index] = tmp_series[group.index].map(normalize_index_phytomer) 
        
    return dimT_dataframe
        
        
    