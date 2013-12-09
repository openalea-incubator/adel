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
Provides functions to construct the :ref:`phenT_abs <phenT_abs>`, the :ref:`phenT_first <phenT_first>`, the :ref:`phenT <phenT>` and the 
:ref:`HS_GL_SSI_T <HS_GL_SSI_T>` dataframes.

The :ref:`phenT_abs <phenT_abs>`, the :ref:`phenT_first <phenT_first>`, the :ref:`phenT <phenT>` and the :ref:`HS_GL_SSI_T <HS_GL_SSI_T>` dataframes are 
described in the User Guide (see :ref:`adel_user`).

Authors: Mariem Abichou, Camille Chambon, Bruno Andrieu
'''

import numpy as np
import pandas

from adel.plantgen import params, tools


def create_phenT_tmp(axeT_tmp, dynT_):
    '''
    Create the *phenT_tmp* dataframe. 
    Compute all the columns, but the column 'TT_del_phytomer' is temporary and will 
    be recalculated in create_phenT_abs using dimT_abs. 
           
    :Parameters:
    
        - `axeT_tmp` (:class:`pandas.DataFrame`) - the :ref:`axeT_tmp <axeT_tmp>` dataframe.
        - `dynT_` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe.
          
    :Returns:
        The *phenT_tmp* dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * *dynT_* must be completely filled, i.e. must not contain any 
          NA value.

    '''
    if not (dynT_.count().max() == dynT_.count().min() == dynT_.index.size):
        raise tools.InputError("dynT contains NA values")
    
    id_phen_list = []
    index_phytomer_list = []
    for (id_phen, N_phytomer_potential), axeT_tmp_group in axeT_tmp.groupby(['id_phen', 'N_phytomer_potential']):
        id_phen_list.extend(np.repeat(id_phen, N_phytomer_potential + 1))
        index_phytomer_list.extend(range(N_phytomer_potential + 1))
        
    phenT_abs = pandas.DataFrame(index=range(len(id_phen_list)), 
                                 columns=['id_phen', 'index_phytomer', 'TT_app_phytomer', 'TT_col_phytomer', 'TT_sen_phytomer', 'TT_del_phytomer'],
                                 dtype=float)
    
    phenT_abs['id_phen'] = id_phen_list
    phenT_abs['index_phytomer'] = index_phytomer_list
    
    phenT_abs_grouped = phenT_abs.groupby('id_phen')
    dynT_grouped = dynT_.groupby(['id_cohort', 'N_phytomer_potential'])
    
    for (id_cohort, N_phytomer_potential, id_phen), axeT_tmp_group in axeT_tmp.groupby(['id_cohort', 'N_phytomer_potential', 'id_phen']):
        phenT_abs_group = phenT_abs_grouped.get_group(id_phen)
        dynT_group = dynT_grouped.get_group((id_cohort, N_phytomer_potential))
        dynT_row = dynT_group.ix[dynT_group[dynT_group['id_axis'] == dynT_group['id_axis'].max()].first_valid_index()]
        a_cohort, TT_col_0, TT_col_break, TT_col_N_phytomer_potential = \
            dynT_row[['a_cohort', 'TT_col_0', 'TT_col_break', 'TT_col_N_phytomer_potential']]
        
        HS_break = a_cohort * (TT_col_break - TT_col_0)
        a2 = (N_phytomer_potential - HS_break) / (TT_col_N_phytomer_potential - TT_col_break)
        
        phenT_abs['TT_col_phytomer'][phenT_abs_group.index] = \
            phenT_abs_group['TT_col_phytomer'] = \
                phenT_abs_group['index_phytomer'].apply(_calculate_TT_col_phytomer, args=(HS_break, TT_col_0, a_cohort, a2, TT_col_break))
        
        first_leaf_indexes = phenT_abs_group.index[0:2]
        phenT_abs['TT_app_phytomer'][first_leaf_indexes] = \
            phenT_abs_group['TT_app_phytomer'][first_leaf_indexes] = \
                phenT_abs_group['TT_col_phytomer'][first_leaf_indexes].apply(
                    _calculate_TT_app_phytomer, args=(TT_col_break, a_cohort, HS_break, N_phytomer_potential, TT_col_N_phytomer_potential, a2, params.DELAIS_PHYLL_COL_TIP_1ST))
        
        other_leaves_indexes = phenT_abs_group.index - phenT_abs_group.index[0:2]
        phenT_abs['TT_app_phytomer'][other_leaves_indexes] = \
            phenT_abs_group['TT_app_phytomer'][other_leaves_indexes] = \
                phenT_abs_group['TT_col_phytomer'][other_leaves_indexes].apply(
                    _calculate_TT_app_phytomer, args=(TT_col_break, a_cohort, HS_break, N_phytomer_potential, TT_col_N_phytomer_potential, a2, params.DELAIS_PHYLL_COL_TIP_NTH))
             
        id_axis, n0, n1, n2, t0, t1, a, c = \
            dynT_row[['id_axis', 'n0', 'n1', 'n2', 't0', 't1', 'a', 'c']]
        HS_1, HS_2, GL_2, GL_3, GL_4 = _calculate_HS_GL_polynomial(HS_break, id_axis, a_cohort, TT_col_0, TT_col_break, TT_col_N_phytomer_potential, n0, n1, n2, t0, t1, a, c, a2)
        
        phenT_abs['TT_sen_phytomer'][phenT_abs_group.index] = \
            phenT_abs_group['TT_sen_phytomer'] = \
                phenT_abs_group['index_phytomer'].apply(_calculate_TT_sen_phytomer, args=(HS_break, HS_1, HS_2, GL_2, GL_3, GL_4, t0, t1, TT_col_N_phytomer_potential, N_phytomer_potential))
        
        phenT_abs['TT_del_phytomer'][phenT_abs_group.index] = \
            phenT_abs_group['TT_del_phytomer'] = \
                _calculate_TT_del_phytomer(a_cohort, phenT_abs_group['TT_sen_phytomer'])
        
    return phenT_abs


def create_phenT_abs(phenT_tmp, axeT_, dimT_abs):
    '''
    Create the :ref:`phenT_abs <phenT_abs>` dataframe.
    
    :Parameters:
    
        - `phenT_tmp` (:class:`pandas.DataFrame`) - the *phenT_tmp* dataframe.
        - `axeT_` (:class:`pandas.DataFrame`) - the :ref:`axeT <axeT>` dataframe.
        - `dimT_abs` (:class:`pandas.DataFrame`) - the :ref:`dimT_abs <dimT_abs>` dataframe.
        
    :Returns:
        The :ref:`phenT_abs <phenT_abs>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
    
    '''
    phenT_abs = phenT_tmp.select(lambda idx: phenT_tmp['id_phen'][idx] in axeT_['id_phen'].values)
    axeT_grouped = axeT_.groupby('id_phen')
    dimT_abs_grouped = dimT_abs.groupby('id_dim')
    for id_phen, phenT_abs_group in phenT_abs.groupby('id_phen'):
        axeT_group = axeT_grouped.get_group(id_phen)
        id_dim = axeT_group['id_dim'].max()
        dimT_abs_group = dimT_abs_grouped.get_group(id_dim)
        non_zero_L_internode_group = dimT_abs_group[dimT_abs_group['L_internode'] != 0]
        min_index_phytomer = non_zero_L_internode_group['index_phytomer'].min()
        indexes_to_ceil = phenT_abs_group[phenT_abs_group['index_phytomer'] >= min_index_phytomer].index
        phenT_abs['TT_del_phytomer'][indexes_to_ceil] = params.TT_DEL_FHAUT
    return phenT_abs


def _calculate_TT_col_phytomer(index_phytomer, HS_break, TT_col_0, a_cohort, a2, TT_col_break):
    if index_phytomer < HS_break: # linear mode
        TT_col_phytomer = TT_col_0 + index_phytomer / a_cohort
    else: # bilinear mode
        TT_col_phytomer = (index_phytomer - HS_break) / a2 + TT_col_break
    return TT_col_phytomer


def _calculate_TT_app_phytomer(TT_col_phytomer, TT_col_break, a_cohort, HS_break, N_phytomer_potential, TT_col_N_phytomer_potential, a2, delais_phyll_col_tip):
    if TT_col_phytomer < TT_col_break:
        TT_app_phytomer = TT_col_phytomer - (delais_phyll_col_tip / a_cohort)
    else:
        TT_app_phytomer = TT_col_phytomer - (delais_phyll_col_tip / a2)
        if TT_app_phytomer <= HS_break:
            TT_app_phytomer = TT_col_break - (delais_phyll_col_tip - a2 * (TT_col_phytomer - TT_col_break)) / a_cohort
    return TT_app_phytomer


def _calculate_TT_sen_phytomer(index_phytomer, HS_break, HS_1, HS_2, GL_2, GL_3, GL_4, t0, t1, TT_col_N_phytomer_potential, N_phytomer_potential):  
    # define HS according to index_phytomer
    if index_phytomer < HS_break: # linear mode
        HS = HS_1
    else: # bilinear mode
        HS = HS_2
    GL_1 = HS
    # Suppose we are in the [0,t0] phase.
    GL = GL_1
    SSI = HS - GL
    if index_phytomer == 0:
        TT_sen_phytomer = t0
    else:
        # Find (SSI - index_phytomer) real root.
        SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
        if SSI_root_array.size == 0 \
            or (SSI_root_array[0] > t0 \
                and not np.allclose(SSI_root_array[0], t0)):
            # Suppose we are in the ]t0,t1] phase.
            GL = GL_2
            SSI = HS - GL
            # Find (SSI - index_phytomer) real root.
            SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
            if SSI_root_array.size == 0 \
                or (SSI_root_array[0] <= t0 \
                    and not np.allclose(SSI_root_array[0], t0)) \
                or (SSI_root_array[0] > t1 \
                    and not np.allclose(SSI_root_array[0], t1)):
                # Suppose we are in the ]t1,TT_col_N_phytomer_potential] phase.
                GL = GL_3
                SSI = HS - GL
                # Find (SSI - index_phytomer) real root.
                SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
                if SSI_root_array.size == 0 \
                    or (SSI_root_array[0] <= t1 \
                        and not np.allclose(SSI_root_array[0], t1)) \
                    or (SSI_root_array[0] > TT_col_N_phytomer_potential \
                        and not np.allclose(SSI_root_array[0], TT_col_N_phytomer_potential)):
                    # We must be in the ]TT_col_N_phytomer_potential,infinity[ phase.
                    GL = GL_4
                    SSI = HS - GL
                    # Find (SSI - index_phytomer) real root.
                    SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
                    if SSI_root_array.size == 0 \
                        or (SSI_root_array[0] <= TT_col_N_phytomer_potential \
                            and not np.allclose(SSI_root_array[0], TT_col_N_phytomer_potential)):
                        raise Exception('ERROR !!!!! This shouldn\'t occurred')
                    if HS(SSI_root_array[0]) > N_phytomer_potential:
                        HS = np.poly1d([N_phytomer_potential])
                    if GL(SSI_root_array[0]) < 0.0:
                        GL = np.poly1d([0.0])
                    SSI = HS - GL
                    # Find (SSI - index_phytomer) real root again.
                    SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
                    if SSI_root_array.size == 0 \
                        or (SSI_root_array[0] <= TT_col_N_phytomer_potential \
                            and not np.allclose(SSI_root_array[0], TT_col_N_phytomer_potential)):
                        raise Exception('ERROR !!!!! This shouldn\'t occurred')
        TT_sen_phytomer = SSI_root_array[0]
    
    return TT_sen_phytomer


def _calculate_HS_GL_polynomial(HS_break, id_axis, a_cohort, TT_col_0, TT_col_break, TT_col_N_phytomer_potential, n0, n1, n2, t0, t1, a, c, a2):
    # define HS(TT)
    HS_1 = np.poly1d([a_cohort, - a_cohort * TT_col_0]) # index_phytomer < HS_break
    HS_2 = np.poly1d([a2, - a2 * TT_col_break + HS_break]) # index_phytomer >= HS_break
    # define GL(TT) for all phases except TT < t0 (because it depends on index_phytomer)
    if id_axis == 'MS':
        GL_2 = np.poly1d([(n1 - n0) / (t1 - t0), n0 - t0 * (n1 - n0) / (t1 - t0)])
        GL_3 = np.poly1d([(n2 - n1) / (TT_col_N_phytomer_potential - t1), n1 - t1 * (n2 - n1) / (TT_col_N_phytomer_potential - t1)])
    else: # tillers
        GL_2 = np.poly1d([n0])
        GL_3 = np.poly1d([(n2 - n0) / (TT_col_N_phytomer_potential - t1), n0 - t1 * (n2 - n0) / (TT_col_N_phytomer_potential - t1)])
    GL_4 = np.poly1d([a, - 3 * a * TT_col_N_phytomer_potential, 3 * a * TT_col_N_phytomer_potential**2 + c, - a * TT_col_N_phytomer_potential**3 - c * TT_col_N_phytomer_potential + n2])
    return HS_1, HS_2, GL_2, GL_3, GL_4


def _calculate_TT_del_phytomer(a_cohort, TT_sen_phytomer_series):
    TT_del_phytomer_series = TT_sen_phytomer_series + params.DELAIS_PHYLL_SEN_DISP / a_cohort
    return TT_del_phytomer_series

    
def create_phenT_first(phenT_abs):
    '''
    Create the :ref:`phenT_first <phenT_first>` dataframe.

    :Parameters:
    
        - `phenT_abs` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
        
    :Returns:
        The :ref:`phenT_first <phenT_first>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * *phenT_abs* must be completely filled, i.e. must not contain 
          any NA value.
    
    '''
    if not (phenT_abs.count().max() == phenT_abs.count().min() == phenT_abs.index.size):
        raise tools.InputError("phenT_abs contains NA values")
    
    phenT_first = phenT_abs.select(lambda idx: True if phenT_abs['index_phytomer'][idx] == 1 else False)
    return phenT_first


def create_phenT(phenT_abs, phenT_first):
    '''
    Create the :ref:`phenT <phenT>` dataframe.

    :Parameters:
    
        - `phenT_abs` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
        - `phenT_first` (:class:`pandas.DataFrame`) - the :ref:`phenT_first <phenT_first>` dataframe.
        
    :Returns:
        The :ref:`phenT <phenT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * *phenT_abs* must be completely filled, i.e. must not contain 
          any NA value.
    
    '''
    if not (phenT_abs.count().max() == phenT_abs.count().min() == phenT_abs.index.size):
        raise tools.InputError("phenT_abs contains NA values")
    
    phenT_ = pandas.DataFrame(index=phenT_abs.index, columns=['id_phen', 'index_rel_phytomer', 'dTT_app_phytomer', 'dTT_col_phytomer', 'dTT_sen_phytomer', 'dTT_del_phytomer'], dtype=float)
    phenT_['id_phen'] = phenT_abs['id_phen']
    tmp_series = pandas.Series(phenT_.index)
    for name, group in phenT_abs.groupby('id_phen'):
        current_first_leaf_row = phenT_first[phenT_first['id_phen'] == name]
        def normalize_index_phytomer(i):
            return group['index_phytomer'][i] / float(group['index_phytomer'][group.index[-1]])
        phenT_['index_rel_phytomer'][group.index] = tmp_series[group.index].map(normalize_index_phytomer)
        def get_relative_TT_col_phytomer(i):
            return group['TT_col_phytomer'][i] - current_first_leaf_row['TT_col_phytomer'][current_first_leaf_row.first_valid_index()]
        phenT_['dTT_col_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_col_phytomer)
        def get_relative_TT_app_phytomer(i):
            return group['TT_app_phytomer'][i] - current_first_leaf_row['TT_app_phytomer'][current_first_leaf_row.first_valid_index()]
        phenT_['dTT_app_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_app_phytomer)
        def get_relative_TT_sen_phytomer(i):
            return group['TT_sen_phytomer'][i] - current_first_leaf_row['TT_sen_phytomer'][current_first_leaf_row.first_valid_index()]
        phenT_['dTT_sen_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_sen_phytomer)
        def get_relative_TT_del_phytomer(i):
            return group['TT_del_phytomer'][i] - current_first_leaf_row['TT_del_phytomer'][current_first_leaf_row.first_valid_index()]
        phenT_['dTT_del_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_del_phytomer)
        
    return phenT_


def create_HS_GL_SSI_T(phenT_abs, axeT_, dynT_):
    '''
    Create the :ref:`HS_GL_SSI_T <HS_GL_SSI_T>` dataframe.

    :Parameters:
    
        - `phenT_abs` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
        - `axeT_` (:class:`pandas.DataFrame`) - the *axeT_* dataframe.
        - `dynT_` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe.
        
    :Returns:
        The :ref:`HS_GL_SSI_T <HS_GL_SSI_T>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    '''
    
    HS_GL_SSI_dynamic = pandas.DataFrame(columns=['id_phen', 'TT', 'HS', 'GL', 'SSI'])
    
    phenT_abs_grouped = phenT_abs.groupby('id_phen')
    dynT_grouped = dynT_.groupby(['id_cohort', 'N_phytomer_potential'])
    for (id_cohort, N_phytomer_potential, id_phen), axeT_group in axeT_.groupby(['id_cohort', 'N_phytomer_potential', 'id_phen']):
        phenT_abs_group = phenT_abs_grouped.get_group(id_phen)
        dynT_group = dynT_grouped.get_group((id_cohort, N_phytomer_potential))
        dynT_row = dynT_group.ix[dynT_group['cardinality'].idxmax()]
        id_axis, a_cohort, TT_col_0, TT_col_break, TT_col_N_phytomer_potential, n0, n1, n2, t0, t1, a, c = \
            dynT_row[['id_axis', 'a_cohort', 'TT_col_0', 'TT_col_break', 'TT_col_N_phytomer_potential', 'n0', 'n1', 'n2', 't0', 't1', 'a', 'c']]
            
        HS_break = a_cohort * (TT_col_break - TT_col_0)
        a2 = (N_phytomer_potential - HS_break) / (TT_col_N_phytomer_potential - TT_col_break)
        
        HS_1, HS_2, GL_2, GL_3, GL_4 = _calculate_HS_GL_polynomial(HS_break, id_axis, a_cohort, TT_col_0, TT_col_break, TT_col_N_phytomer_potential, n0, n1, n2, t0, t1, a, c, a2)
        
        t0, t1, TT_col_N_phytomer_potential, TT_col_break = np.round([t0, t1, TT_col_N_phytomer_potential, TT_col_break]).astype(int)
        
        TT_1_1 = np.arange(0, t0)
        TT_1_2 = np.arange(t0, t0)
        TT_2_1 = np.arange(t0, t1)
        TT_2_2 = np.arange(t1, t1)
        TT_3_1 = np.arange(t1, TT_col_N_phytomer_potential)
        TT_3_2 = np.arange(TT_col_N_phytomer_potential, TT_col_N_phytomer_potential)
        TT_4_1 = np.arange(TT_col_N_phytomer_potential, params.TT_DEL_FHAUT)
        TT_4_2 = np.arange(params.TT_DEL_FHAUT, params.TT_DEL_FHAUT)
        
        if TT_col_break != 0.0: # bilinear mode
            if TT_col_break <= t0:
                TT_1_1 = np.arange(0, TT_col_break)
                TT_1_2 = np.arange(TT_col_break, t0)
            elif TT_col_break <= t1:
                TT_2_1 = np.arange(t0, TT_col_break)
                TT_2_2 = np.arange(TT_col_break, t1)
            elif TT_col_break <= TT_col_N_phytomer_potential:
                TT_3_1 = np.arange(t1, TT_col_break)
                TT_3_2 = np.arange(TT_col_break, TT_col_N_phytomer_potential)
            else:
                TT_4_1 = np.arange(TT_col_N_phytomer_potential, TT_col_break)
                TT_4_2 = np.arange(TT_col_break, params.TT_DEL_FHAUT)
            
        HS_1_TT_1_1 = np.clip(HS_1(TT_1_1), 0.0, N_phytomer_potential)
        HS_1_TT_2_1 = np.clip(HS_1(TT_2_1), 0.0, N_phytomer_potential)
        HS_1_TT_3_1 = np.clip(HS_1(TT_3_1), 0.0, N_phytomer_potential)
        HS_1_TT_4_1 = np.clip(HS_1(TT_4_1), 0.0, N_phytomer_potential)
        HS_1_TT_1_2 = np.clip(HS_1(TT_1_2), 0.0, N_phytomer_potential)
        HS_1_TT_2_2 = np.clip(HS_1(TT_2_2), 0.0, N_phytomer_potential)
        HS_1_TT_3_2 = np.clip(HS_1(TT_3_2), 0.0, N_phytomer_potential)
        HS_1_TT_4_2 = np.clip(HS_1(TT_4_2), 0.0, N_phytomer_potential)
        
        GL_1_TT_1_1 = np.clip(HS_1(TT_1_1), 0.0, 1000.0)
        GL_2_TT_2_1 = np.clip(GL_2(TT_2_1), 0.0, 1000.0)
        GL_3_TT_3_1 = np.clip(GL_3(TT_3_1), 0.0, 1000.0)
        GL_4_TT_4_1 = np.clip(GL_4(TT_4_1), 0.0, 1000.0)
        GL_1_TT_1_2 = np.clip(HS_2(TT_1_2), 0.0, 1000.0)
        GL_2_TT_2_2 = np.clip(GL_2(TT_2_2), 0.0, 1000.0)
        GL_3_TT_3_2 = np.clip(GL_3(TT_3_2), 0.0, 1000.0)
        GL_4_TT_4_2 = np.clip(GL_4(TT_4_2), 0.0, 1000.0)
        
        SSI_1_TT_1_1 = HS_1_TT_1_1 - GL_1_TT_1_1
        SSI_2_TT_2_1 = HS_1_TT_2_1 - GL_2_TT_2_1
        SSI_3_TT_3_1 = HS_1_TT_3_1 - GL_3_TT_3_1
        SSI_4_TT_4_1 = HS_1_TT_4_1 - GL_4_TT_4_1
        SSI_1_TT_1_2 = HS_1_TT_1_2 - GL_1_TT_1_2
        SSI_2_TT_2_2 = HS_1_TT_2_2 - GL_2_TT_2_2
        SSI_3_TT_3_2 = HS_1_TT_3_2 - GL_3_TT_3_2
        SSI_4_TT_4_2 = HS_1_TT_4_2 - GL_4_TT_4_2
        
        HS_GL_SSI_dynamic_group = pandas.DataFrame(index=np.arange(params.TT_DEL_FHAUT), columns=HS_GL_SSI_dynamic.columns)
        HS_GL_SSI_dynamic_group['id_phen'] = np.repeat(id_phen, HS_GL_SSI_dynamic_group.index.size)
        HS_GL_SSI_dynamic_group['TT'] = pandas.Series(np.concatenate((TT_1_1, TT_1_2, TT_2_1, TT_2_2, TT_3_1, TT_3_2, TT_4_1, TT_4_2)))
        HS_GL_SSI_dynamic_group['HS'] = pandas.Series(np.concatenate((HS_1_TT_1_1, HS_1_TT_1_2, HS_1_TT_2_1, HS_1_TT_2_2, HS_1_TT_3_1, HS_1_TT_3_2, HS_1_TT_4_1, HS_1_TT_4_2)))
        HS_GL_SSI_dynamic_group['GL'] = pandas.Series(np.concatenate((GL_1_TT_1_1, GL_1_TT_1_2, GL_2_TT_2_1, GL_2_TT_2_2, GL_3_TT_3_1, GL_3_TT_3_2, GL_4_TT_4_1, GL_4_TT_4_2)))
        HS_GL_SSI_dynamic_group['SSI'] = pandas.Series(np.concatenate((SSI_1_TT_1_1, SSI_1_TT_1_2, SSI_2_TT_2_1, SSI_2_TT_2_2, SSI_3_TT_3_1, SSI_3_TT_3_2, SSI_4_TT_4_1, SSI_4_TT_4_2)))
        
        HS_GL_SSI_dynamic = HS_GL_SSI_dynamic.append(HS_GL_SSI_dynamic_group, ignore_index=True)
    
    return HS_GL_SSI_dynamic


    