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

from adel.plantgen import params


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
    id_dim_tmp_list = _gen_id_dim_tmp_list(dynT_tmp)
    index_phytomer_tmp_list = _gen_index_phytomer_tmp_list(id_dim_tmp_list)
    L_blade_list = [np.nan for i in range(len(id_dim_tmp_list))]
    W_blade_list = [np.nan for i in range(len(id_dim_tmp_list))]
    L_sheath_list = [np.nan for i in range(len(id_dim_tmp_list))]
    W_sheath_list = [np.nan for i in range(len(id_dim_tmp_list))]
    L_internode_list = [np.nan for i in range(len(id_dim_tmp_list))]
    W_internode_list = [np.nan for i in range(len(id_dim_tmp_list))]
    dimT_tmp_array = np.array([id_dim_tmp_list, index_phytomer_tmp_list, L_blade_list, W_blade_list, L_sheath_list, W_sheath_list, L_internode_list, W_internode_list]).transpose()
    return pandas.DataFrame(dimT_tmp_array, columns=['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode'], dtype=float)
    

def create_dimT_abs(axeT, dimT_user, phenT_abs):
    '''
    Create the :ref:`dimT_abs <dimT_abs>` dataframe filling the *dimT_user* dataframe.

    :Parameters:
        
        - `axeT` (:class:`pandas.DataFrame`) - the :ref:`axeT <axeT>` dataframe.
        - `dimT_user` (:class:`pandas.DataFrame`) - the :ref:`dimT <dimT>` dataframe set by the user.
        - `phenT_abs` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
    
    :Returns: 
        The :ref:`dimT_abs <dimT_abs>` dataframe.
    
    :Returns Type: 
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * *axeT*, *dimT_user* and *phenT_abs* must be of type :class:`pandas.DataFrame`
        * in *dimT_user*, the column *id_dim* and *index_phytomer* must be 
          completely filled, i.e. they must not contain any NA value.
        * in *dimT_user*, the rows which describe the first axis must be 
          completely filled, i.e. there must not contain any NA value.
        
    '''
    assert isinstance(axeT, pandas.DataFrame) and \
           isinstance(dimT_user, pandas.DataFrame) and \
           isinstance(phenT_abs, pandas.DataFrame)
    assert dimT_user['id_dim'].count() == dimT_user['id_dim'].size
    assert dimT_user['index_phytomer'].count() == dimT_user['index_phytomer'].size
    first_axis_rows_number = int(str(int(dimT_user['id_dim'][0]))[-2:])
    for column_name in dimT_user:
        assert dimT_user[column_name][:first_axis_rows_number].map(lambda x: x == 0.0 and 1.0 or x).fillna(0.0).all()
           
    dimT_abs_dataframe = _init_dimT_abs(axeT['id_dim'].unique(), 
                                        dimT_user, 
                                        phenT_abs)
    
    L_blade_is_null = dimT_abs_dataframe['L_blade'].isnull()
    row_indexes_to_fit = L_blade_is_null[L_blade_is_null == True].index
    
    _gen_lengths(first_axis_rows_number, row_indexes_to_fit, dimT_abs_dataframe)
    
    _gen_widths(first_axis_rows_number, row_indexes_to_fit, dimT_abs_dataframe)
            
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


def _init_dimT_abs(id_dim_from_axeT, dimT_user, phenT_abs):
    '''Initialize dimT_abs.'''
    dimT_user_copy = dimT_user.copy()
    dimT_user_copy['TT_em_phytomer'] = pandas.Series(phenT_abs.drop(phenT_abs.groupby('index_phytomer').groups[0.0])['TT_em_phytomer'].values)
    dimT_user_copy['id_dim_old'] = dimT_user_copy['id_dim'].astype(int)
    dimT_user_copy['id_dim'] = dimT_user_copy['id_dim_old'] * 10 + 1
    dimT_user_copy['is_ear'] = 1
    idx = dimT_user_copy.index.tolist()
    idx.sort(reverse=True)
    dimT_user_copy['idx'] = idx
    id_dim_from_axeT = id_dim_from_axeT.astype(int)
    id_dim_from_axeT_str = np.char.mod('%d', id_dim_from_axeT)
    id_dim_from_axeT_3int = np.core.defchararray.ljust(id_dim_from_axeT_str, 3).astype(int)
    dimT_abs = pandas.DataFrame()
    for id_dim_from_dimT_user, dimT_user_group in dimT_user_copy.groupby('id_dim_old'):
        if id_dim_from_dimT_user in id_dim_from_axeT_3int:
            dimT_abs = pandas.concat([dimT_abs, dimT_user_group], ignore_index=True)
            id_dim_regress = id_dim_from_dimT_user * 10
            if id_dim_regress in id_dim_from_axeT:
                regress_dimT_group = dimT_user_group.copy()
                regress_dimT_group['id_dim'] = id_dim_regress
                regress_dimT_group[['L_blade', 'L_internode', 'L_sheath', 'W_blade', 'W_internode', 'W_sheath']] = np.nan
                regress_dimT_group['is_ear'] = 0
                dimT_abs = pandas.concat([dimT_abs, regress_dimT_group], ignore_index=True)
    
    dimT_abs = dimT_abs.sort_index(by=['is_ear', 'idx'], ascending=False)
    dimT_abs = dimT_abs.drop(['idx', 'id_dim_old'], axis=1)
    
    dimT_abs_array = np.array([dimT_abs['id_dim'].values, 
                               dimT_abs['index_phytomer'].values,
                               dimT_abs['L_blade'].values,
                               dimT_abs['W_blade'].values,
                               dimT_abs['L_sheath'].values,
                               dimT_abs['W_sheath'].values,
                               dimT_abs['L_internode'].values,
                               dimT_abs['W_internode'].values,
                               dimT_abs['TT_em_phytomer'].values,
                               dimT_abs['is_ear'].values]).transpose()
    dimT_abs = pandas.DataFrame(dimT_abs_array, 
                                columns=['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode', 'TT_em_phytomer', 'is_ear'])

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


def _gen_lengths(first_axis_rows_number, row_indexes_to_fit, dimT_abs_dataframe):
    '''Fit the lengths in-place.'''
    TT_em_phytomer_series = dimT_abs_dataframe['TT_em_phytomer']
    MS_last_TT_em_phytomer = TT_em_phytomer_series[first_axis_rows_number - 1]
    indexes_to_ceil = TT_em_phytomer_series[TT_em_phytomer_series > MS_last_TT_em_phytomer].index
    indexes_to_ceil = indexes_to_ceil.intersection(row_indexes_to_fit)
    
    for length in ['L_blade', 'L_sheath', 'L_internode']:
        first_axis_rows_series = dimT_abs_dataframe[length][:first_axis_rows_number]
        first_null_axis_rows_series = first_axis_rows_series[first_axis_rows_series == 0.0]
        position_of_first_non_null_data = first_null_axis_rows_series.index.size
        polynomial_coefficients_array = np.polyfit(dimT_abs_dataframe['TT_em_phytomer'][position_of_first_non_null_data:first_axis_rows_number].values, 
                                                   dimT_abs_dataframe[length][position_of_first_non_null_data:first_axis_rows_number].values, 6)
        for name, group in dimT_abs_dataframe.ix[row_indexes_to_fit].groupby(by='id_dim'):
            TT_em_phytomer_group = group['TT_em_phytomer']
            fitted_length_series = dimT_abs_dataframe[length]
            fitted_length_series[group.index] = np.polyval(polynomial_coefficients_array, 
                                                           TT_em_phytomer_group)
            if group['is_ear'][group.first_valid_index()] == 0:
                # regression
                for i in range(len(params.REGRESSION_OF_DIMENSIONS[length])):
                    fitted_length_series[group.index[i-3]] *= (1.0 - params.REGRESSION_OF_DIMENSIONS[length][i])

        # thresholding
        if length == 'L_internode':
            MS_first_TT_em_phytomer = TT_em_phytomer_series[position_of_first_non_null_data]
            indexes_to_threshold = TT_em_phytomer_series[TT_em_phytomer_series <= MS_first_TT_em_phytomer].index
            indexes_to_threshold = indexes_to_threshold.intersection(row_indexes_to_fit)
            dimT_abs_dataframe[length][indexes_to_threshold] = 0.0
        # ceiling
        MS_last_length = dimT_abs_dataframe[length][first_axis_rows_number - 1]
        dimT_abs_dataframe[length][indexes_to_ceil] = MS_last_length


def _gen_widths(first_axis_rows_number, row_indexes_to_fit, dimT_abs_dataframe):
    '''Fit the widths in-place.'''
    for width in ['W_blade', 'W_sheath', 'W_internode']:
        first_axis_rows_series = dimT_abs_dataframe[width][:first_axis_rows_number]
        first_null_axis_rows_series = first_axis_rows_series[first_axis_rows_series == 0.0]
        position_of_first_non_null_data = first_null_axis_rows_series.index.size
        first_axis_W_tuple = (dimT_abs_dataframe[width][position_of_first_non_null_data], dimT_abs_dataframe[width][first_axis_rows_number - 1])
        if width == 'W_internode':
            first_axis_TT_em_phytomer_tuple = (dimT_abs_dataframe['TT_em_phytomer'][position_of_first_non_null_data], dimT_abs_dataframe['TT_em_phytomer'][first_axis_rows_number - 1])
        else:
            first_axis_TT_em_phytomer_tuple = None
        for name, group in dimT_abs_dataframe.ix[row_indexes_to_fit].groupby(by='id_dim'):
            TT_em_phytomer_sub_series = group['TT_em_phytomer']
            fitted_width_series = dimT_abs_dataframe[width]
            if first_axis_TT_em_phytomer_tuple is None: # NOT internode
                x1 = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series.first_valid_index()]
                x2 = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series.last_valid_index()]
            else: # internode
                # get TT_em_phytomer of the main stem first phytomer which has a 
                # TT_em_phytomer greater than the first main stem phytomer with 
                # a non null width
                valid_TT_em_phytomers = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series >= first_axis_TT_em_phytomer_tuple[0]]
                x1 = valid_TT_em_phytomers[valid_TT_em_phytomers.index[0]]
                # get TT_em_phytomer of the main stem last phytomer which has a 
                # TT_em_phytomer lesser than the last main stem phytomer
                valid_TT_em_phytomers = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series <= first_axis_TT_em_phytomer_tuple[-1]]
                x2 = valid_TT_em_phytomers[valid_TT_em_phytomers.index[-1]]
                
            y1 = first_axis_W_tuple[0]
            y2 = first_axis_W_tuple[1]
            polynomial_coefficient_array = np.polyfit(np.array([x1, x2]), np.array([y1, y2]), 1)
            fitted_width_series[group.index] = np.polyval(polynomial_coefficient_array, TT_em_phytomer_sub_series)
            if group['is_ear'][group.first_valid_index()] == 0:
                # regression
                for i in range(len(params.REGRESSION_OF_DIMENSIONS[width])):
                    fitted_width_series[group.index[i-3]] *= (1.0 - params.REGRESSION_OF_DIMENSIONS[width][i])
        
        if width == 'W_internode':
            # thresholding of the width of the phytomers which have a TT_em_phytomer 
            # lesser than the TT_em_phytomer of the first main stem which has a non null width
            MS_first_TT_em_phytomer = first_axis_TT_em_phytomer_tuple[0]
            TT_em_phytomer_sub_series = dimT_abs_dataframe['TT_em_phytomer']
            indexes_to_threshold = TT_em_phytomer_sub_series[TT_em_phytomer_sub_series <= MS_first_TT_em_phytomer].index
            indexes_to_threshold = indexes_to_threshold.intersection(row_indexes_to_fit)
            dimT_abs_dataframe[width][indexes_to_threshold] = 0.0
            # ceiling of the the width of the phytomers which have a TT_em_phytomer 
            # greater than the TT_em_phytomer of the last main stem 
            MS_last_TT_em_phytomer = first_axis_TT_em_phytomer_tuple[1]
            width_series = dimT_abs_dataframe[width]
            indexes_to_ceil = width_series[width_series > MS_last_TT_em_phytomer].index
            indexes_to_ceil = indexes_to_ceil.intersection(row_indexes_to_fit)
            dimT_abs_dataframe[width][indexes_to_ceil] = first_axis_W_tuple[1]
        dimT_abs_dataframe[width] = dimT_abs_dataframe[width].clip_lower(0.0)
    

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
        
        
    