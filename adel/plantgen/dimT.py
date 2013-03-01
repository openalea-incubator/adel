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
    
    .. warning:: *dynT_tmp* must be a :class:`pandas.DataFrame`.
    
    '''
    assert isinstance(dynT_tmp, pandas.DataFrame)
    id_dim_list = _gen_id_dim_list(dynT_tmp)
    index_phytomer_list = _gen_index_phytomer_list(id_dim_list)
    L_blade_list = [np.nan for i in range(len(id_dim_list))]
    W_blade_list = [np.nan for i in range(len(id_dim_list))]
    L_sheath_list = [np.nan for i in range(len(id_dim_list))]
    W_sheath_list = [np.nan for i in range(len(id_dim_list))]
    L_internode_list = [np.nan for i in range(len(id_dim_list))]
    W_internode_list = [np.nan for i in range(len(id_dim_list))]
    dimT_tmp_array = np.array([id_dim_list, index_phytomer_list, L_blade_list, W_blade_list, L_sheath_list, W_sheath_list, L_internode_list, W_internode_list]).transpose()
    return pandas.DataFrame(dimT_tmp_array, columns=['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'], dtype=float)
    

def create_dimT_abs(dimT_user, phenT_abs):
    '''
    Create the :ref:`dimT_abs <dimT_abs>` dataframe filling the *dimT_user* dataframe.

    :Parameters:
    
        - `dimT_user` (:class:`pandas.DataFrame`) - the :ref:`dimT <dimT>` dataframe set by the user.
        - `phenT_abs` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
    
    :Returns: 
        The :ref:`dimT_abs <dimT_abs>` dataframe.
    
    :Returns Type: 
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * *dimT_user* and *phenT_abs* must be of type :class:`pandas.DataFrame`
        * in *dimT_user*, the column *id_dim* and *index_phytomer* must be 
          completely filled, i.e. they must not contain any NA value.
        * in *dimT_user*, the rows which describe the first axis must be 
          completely filled, i.e. there must not contain any NA value.
        
    '''
    assert isinstance(dimT_user, pandas.DataFrame) and isinstance(phenT_abs, pandas.DataFrame)
    assert dimT_user['id_dim'].count() == dimT_user['id_dim'].size
    assert dimT_user['index_phytomer'].count() == dimT_user['index_phytomer'].size
    first_axis_rows_number = int(str(int(dimT_user['id_dim'][0]))[-2:])
    for column_name in dimT_user:
        assert dimT_user[column_name][:first_axis_rows_number].map(lambda x: x == 0.0 and 1.0 or x).fillna(0.0).all()
    first_row_number_to_fit = dimT_user['L_blade'].last_valid_index() + 1
    dimT_abs_dataframe = dimT_user.copy()
    dimT_abs_dataframe['TT_em_phytomer'] = pandas.Series(phenT_abs.drop(phenT_abs.groupby('index_phytomer').groups[0.0])['TT_em_phytomer'].values)
    dimT_abs_dataframe['L_blade'] = _gen_length(dimT_abs_dataframe['TT_em_phytomer'], dimT_abs_dataframe['L_blade'], first_axis_rows_number, first_row_number_to_fit)
    dimT_abs_dataframe['L_sheath'] = _gen_length(dimT_abs_dataframe['TT_em_phytomer'], dimT_abs_dataframe['L_sheath'], first_axis_rows_number, first_row_number_to_fit)
    dimT_abs_dataframe['L_internode'] = _gen_length(dimT_abs_dataframe['TT_em_phytomer'], dimT_abs_dataframe['L_internode'], first_axis_rows_number, first_row_number_to_fit, True)
    
    W_internode_first_axis_rows_series = dimT_abs_dataframe['W_internode'][:first_axis_rows_number]
    first_null_axis_rows_series = W_internode_first_axis_rows_series[W_internode_first_axis_rows_series == 0.0]
    position_of_first_non_null_data = first_null_axis_rows_series.index.size
    first_axis_W_blade_tuple = (dimT_abs_dataframe['W_blade'][0], dimT_abs_dataframe['W_blade'][first_axis_rows_number - 1])
    first_axis_W_sheath_tuple = (dimT_abs_dataframe['W_sheath'][0], dimT_abs_dataframe['W_sheath'][first_axis_rows_number - 1])
    first_axis_W_internode_tuple = (dimT_abs_dataframe['W_internode'][position_of_first_non_null_data], dimT_abs_dataframe['W_internode'][first_axis_rows_number - 1])
    first_axis_TT_em_phytomer_tuple = (dimT_abs_dataframe['TT_em_phytomer'][position_of_first_non_null_data], dimT_abs_dataframe['TT_em_phytomer'][first_axis_rows_number - 1])
    for name, group in dimT_abs_dataframe[first_row_number_to_fit:].groupby(by='id_dim'):
        dimT_abs_dataframe['W_blade'][group.index] = _gen_width(group['TT_em_phytomer'], first_axis_TT_em_phytomer_tuple, first_axis_W_blade_tuple)
        dimT_abs_dataframe['W_sheath'][group.index] = _gen_width(group['TT_em_phytomer'], first_axis_TT_em_phytomer_tuple, first_axis_W_sheath_tuple)
        dimT_abs_dataframe['W_internode'][group.index] = _gen_width(group['TT_em_phytomer'], first_axis_TT_em_phytomer_tuple, first_axis_W_internode_tuple, True)
    
    return dimT_abs_dataframe.drop(['TT_em_phytomer'], axis=1)


def _gen_id_dim_list(dyn_tmp_dataframe):
    '''Generate the *id_dim* column.'''
    sorted_id_axis = dyn_tmp_dataframe['id_axis']
    id_dim_list = []
    for id_dim in sorted_id_axis:
        N_phyt = int(str(int(id_dim))[-2:])
        for i in range(N_phyt):
            id_dim_list.append(id_dim)
    return id_dim_list


def _gen_index_phytomer_list(id_dim_list):
    '''Generate the *index_phytomer* column.'''
    index_phytomer_list = []
    i = 0
    while i < len(id_dim_list):
        N_phyt = int(str(int(id_dim_list[i]))[-2:])
        for j in range(1, N_phyt + 1):
            index_phytomer_list.append(j)    
        i = i + j
    return index_phytomer_list


def _gen_length(TT_em_phytomer_series, L_series, first_axis_rows_number, first_row_number_to_fit, is_internode=False):
    '''Generate the *L_\** columns.'''
    fitted_length_series = L_series.copy()
    if is_internode:
        first_axis_rows_series = L_series[:first_axis_rows_number]
        first_null_axis_rows_series = first_axis_rows_series[first_axis_rows_series == 0.0]
        position_of_first_non_null_data = first_null_axis_rows_series.index.size
        polynomial_coefficient_array = np.polyfit(TT_em_phytomer_series[position_of_first_non_null_data:first_axis_rows_number].values, L_series[position_of_first_non_null_data:first_axis_rows_number].values, 6)
        fitted_length_series[first_row_number_to_fit:] = np.polyval(polynomial_coefficient_array, TT_em_phytomer_series[first_row_number_to_fit:])
        MS_first_TT_em_phytomer = TT_em_phytomer_series[position_of_first_non_null_data]
        indexes_to_change = TT_em_phytomer_series[TT_em_phytomer_series <= MS_first_TT_em_phytomer].index
        indexes_to_change = indexes_to_change[indexes_to_change >= first_row_number_to_fit]
        fitted_length_series[indexes_to_change] = 0.0
    else:
        polynomial_coefficient_array = np.polyfit(TT_em_phytomer_series[:first_axis_rows_number].values, L_series[:first_axis_rows_number].values, 6)
        fitted_length_series[first_row_number_to_fit:] = np.polyval(polynomial_coefficient_array, TT_em_phytomer_series[first_row_number_to_fit:])
    MS_last_TT_em_phytomer = TT_em_phytomer_series[first_axis_rows_number - 1]
    indexes_to_change = TT_em_phytomer_series[TT_em_phytomer_series > MS_last_TT_em_phytomer].index
    indexes_to_change = indexes_to_change[indexes_to_change >= first_row_number_to_fit]
    fitted_length_series[indexes_to_change] = L_series[first_axis_rows_number - 1]
    return fitted_length_series


def _gen_width(TT_em_phytomer_sub_series, first_axis_TT_em_phytomer_tuple, first_axis_W_tuple, is_internode=False):
    '''Generate the *W_\** columns.'''
    if is_internode:
        # TT_em_phytomer of the main stem first phytomer which has a TT_em_phytomer greater than the first main stem phytomer with a non null width
        valid_TT_em_phytomers = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series >= first_axis_TT_em_phytomer_tuple[0]]
        x1 = valid_TT_em_phytomers[valid_TT_em_phytomers.index[0]]
        # TT_em_phytomer of the main stem last phytomer which has a TT_em_phytomer lesser than the last main stem phytomer
        valid_TT_em_phytomers = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series <= first_axis_TT_em_phytomer_tuple[-1]]
        x2 = valid_TT_em_phytomers[valid_TT_em_phytomers.index[-1]]
    else:
        x1 = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series.index[0]]
        x2 = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series.index[-1]]
    y1 = first_axis_W_tuple[0]
    y2 = first_axis_W_tuple[1]
    polynomial_coefficient_array = np.polyfit(np.array([x1, x2]), np.array([y1, y2]), 1)
    fitted_width_series = pandas.Series(np.polyval(polynomial_coefficient_array, TT_em_phytomer_sub_series), index=TT_em_phytomer_sub_series.index)
    if is_internode:
        # set to 0.0 the width of the phytomers which have a TT_em_phytomer lesser than the TT_em_phytomer of the first main stem which has a non null width
        MS_first_TT_em_phytomer = first_axis_TT_em_phytomer_tuple[0]
        indexes_to_change = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series <= MS_first_TT_em_phytomer].index
        fitted_width_series[indexes_to_change] = 0.0
        # set to first_axis_TT_em_phytomer_tuple[1] the width of the phytomers which have a TT_em_phytomer greater than the TT_em_phytomer of the last main stem
        MS_last_TT_em_phytomer = first_axis_TT_em_phytomer_tuple[1]
        indexes_to_change = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series > MS_last_TT_em_phytomer].index
        fitted_width_series[indexes_to_change] = first_axis_W_tuple[1]
    return fitted_width_series.clip_lower(0.0)


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
        
        * *dimT_abs* must be a :class:`pandas.DataFrame`.
        * *dimT_abs* must be completely filled, i.e. must not contain any 
          NA value.
                 
    '''
    assert isinstance(dimT_abs, pandas.DataFrame)
    assert dimT_abs.count().max() == dimT_abs.count().min() == dimT_abs.index.size
    dimT_dataframe = pandas.DataFrame(index=dimT_abs.index, columns=['id_dim', 'index_rel_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'])
    dimT_dataframe[['id_dim', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']] = dimT_abs[['id_dim', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']]
    tmp_series = pandas.Series(dimT_dataframe.index)
    for name, group in dimT_abs.groupby('id_dim'):
        def normalize_index_phytomer(i):
            return group['index_phytomer'][i] / float(group['index_phytomer'][group.index[-1]])
        dimT_dataframe['index_rel_phytomer'][group.index] = tmp_series[group.index].map(normalize_index_phytomer) 
        
    return dimT_dataframe
        
        
    