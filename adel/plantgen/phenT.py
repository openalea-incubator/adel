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


def create_phenT_abs(dynT_dataframe, decimal_elongated_internode_number):
    '''
    Create the :ref:`phenT_abs <phenT_abs>` dataframe.

    :Parameters:
    
        - `dynT_dataframe` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe.
        - `decimal_elongated_internode_number` (:class:`float`) - the number of elongated 
          internodes.
        
    :Returns:
        The :ref:`phenT_abs <phenT_abs>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
    
    .. warning:: 
    
        * *dynT_dataframe* must be completely filled, i.e. must not contain any 
          NA value.
    
    '''
    
    tools.checkValidity(dynT_dataframe.count().max() == dynT_dataframe.count().min() == dynT_dataframe.index.size)
    id_phen_list = _gen_id_phen_list(dynT_dataframe)
    absolute_index_phytomer_list = _gen_absolute_index_phytomer_list(id_phen_list)
    absolute_TT_col_phytomer_list = _gen_absolute_TT_col_phytomer_list(dynT_dataframe)
    absolute_TT_em_phytomer_list = _gen_absolute_TT_em_phytomer_list(absolute_TT_col_phytomer_list, dynT_dataframe)
    absolute_TT_sen_phytomer_list = _gen_absolute_TT_sen_phytomer_list(dynT_dataframe)
    absolute_TT_del_phytomer_list = _gen_absolute_TT_del_phytomer_list(id_phen_list, absolute_TT_sen_phytomer_list, dynT_dataframe, decimal_elongated_internode_number)
    phenT_abs_array = np.array([id_phen_list, absolute_index_phytomer_list, absolute_TT_em_phytomer_list, absolute_TT_col_phytomer_list, absolute_TT_sen_phytomer_list, absolute_TT_del_phytomer_list]).transpose()
    return pandas.DataFrame(phenT_abs_array, columns=['id_phen', 'index_phytomer', 'TT_em_phytomer', 'TT_col_phytomer', 'TT_sen_phytomer', 'TT_del_phytomer'])


def _gen_id_phen_list(dynT_dataframe):
    '''Generate the *id_phen* column.'''
    sorted_id_axis = dynT_dataframe['id_axis']
    id_phen_list = []
    for id_phen in sorted_id_axis:
        N_phyt = int(str(int(id_phen))[-2:])
        for i in range(N_phyt + 1):
            id_phen_list.append(id_phen)
    return id_phen_list


def _gen_absolute_index_phytomer_list(id_phen_list):
    '''Generate the *index_phytomer* column.'''
    absolute_index_phytomer_list = []
    i = 0
    while i < len(id_phen_list):
        N_phyt = int(str(int(id_phen_list[i]))[-2:])
        for j in range(N_phyt + 1):
            absolute_index_phytomer_list.append(j)    
        i = i + j + 1
    return absolute_index_phytomer_list


def _gen_absolute_TT_col_phytomer_list(dynT_dataframe):
    '''Generate the *TT_col_phytomer* column.'''
    TT_col_phytomer_list = []
    for i in dynT_dataframe.index:
        N_cohort_i, id_axis_i, cardinality_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, dTT_MS_cohort_i, n0_i, n1_i, n2_i, t0_i, t1_i, hs_t1_i, a_i, c_i, RMSE_gl = dynT_dataframe.ix[i].tolist()
        phytomer_indexes = range(int(Nff_i + 1))
        HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
        a2_i = (Nff_i - HS_break_i) / (TT_col_nff_i - TT_col_break_i)
        for j in phytomer_indexes: 
            if j < HS_break_i: # linear mode
                TT_col_phytomer_i_j = TT_col_0_i + j / a_cohort_i
            else: # bilinear mode
                TT_col_phytomer_i_j = (j - HS_break_i) / a2_i + TT_col_break_i
            TT_col_phytomer_list.append(TT_col_phytomer_i_j) 
    
    return TT_col_phytomer_list
    
    
def _gen_absolute_TT_em_phytomer_list(TT_col_phytomer_list, dynT_dataframe):
    '''Generate the *TT_em_phytomer* column.'''
    TT_em_phytomer_list = []
    current_TT_em_phytomer_row_index = 0
    for i in dynT_dataframe.index:
        N_cohort_i, id_axis_i, cardinality_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, dTT_MS_cohort_i, n0_i, n1_i, n2_i, t0_i, t1_i, hs_t1_i, a_i, c_i, RMSE_gl = dynT_dataframe.ix[i].tolist()
        phytomer_indexes = np.arange(int(Nff_i) + 1) + current_TT_em_phytomer_row_index
        j = 0
        for j in phytomer_indexes: 
            TT_col_phytomer_j = TT_col_phytomer_list[j]
            if TT_col_phytomer_j < TT_col_break_i:
                TT_em_phytomer_j = TT_col_phytomer_j - (params.DELAIS_PHYLL_COL_TIP / a_cohort_i)
            else:
                HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
                a2_i = (Nff_i - HS_break_i) / (TT_col_nff_i - TT_col_break_i)
                tmp_res = TT_col_phytomer_j - (params.DELAIS_PHYLL_COL_TIP / a2_i)
                if tmp_res > HS_break_i:
                    TT_em_phytomer_j = tmp_res
                else:
                    TT_em_phytomer_j = TT_col_break_i - (params.DELAIS_PHYLL_COL_TIP - a2_i * (TT_col_phytomer_j - TT_col_break_i)) / a_cohort_i
            TT_em_phytomer_list.append(TT_em_phytomer_j)
        current_TT_em_phytomer_row_index += phytomer_indexes.size
            
    return TT_em_phytomer_list


def _gen_absolute_TT_sen_phytomer_list(dynT_dataframe):
    '''Generate the *TT_sen_phytomer* column.'''
    TT_sen_phytomer_list = []
    
    def get_real_roots(poly):
        roots_array = poly.r
        null_imag_solutions = filter(lambda x: x.imag == 0.0, roots_array)
        real_solutions = map(lambda x: x.real, null_imag_solutions)
        return np.array(real_solutions)
        
    for i in dynT_dataframe.index:
        TT_sen_phytomer_i_list = []        
        N_cohort_i, id_axis_i, cardinality_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, dTT_MS_cohort_i, n0_i, n1_i, n2_i, t0_i, t1_i, hs_t1_i, a_i, c_i, RMSE_gl = dynT_dataframe.ix[i].tolist()
        HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
        HS_1, HS_2, GL_2, GL_3, GL_4 = _gen_HS_GL_polynomial(HS_break_i, N_cohort_i, id_axis_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, n0_i, n1_i, n2_i, t0_i, t1_i, a_i, c_i)
        phytomer_indexes_i = np.arange(int(Nff_i) + 1)

        for j in phytomer_indexes_i:
            # define HS according to j
            if j < HS_break_i: # linear mode
                HS = HS_1
            else: # bilinear mode
                HS = HS_2
            GL_1 = HS
            # Suppose we are in the [0,t0_i] phase.
            GL = GL_1
            SSI = HS - GL
            if j == 0:
                TT_sen_phytomer_i_list.append(t0_i)
            else:
                # Find (SSI - j) real root.
                SSI_root_array = get_real_roots(SSI - j)
                if SSI_root_array.size == 0 \
                    or (SSI_root_array[0] > t0_i \
                        and not np.allclose(SSI_root_array[0], t0_i)):
                    # Suppose we are in the ]t0_i,t1_i] phase.
                    GL = GL_2
                    SSI = HS - GL
                    # Find (SSI - j) real root.
                    SSI_root_array = get_real_roots(SSI - j)
                    if SSI_root_array.size == 0 \
                        or (SSI_root_array[0] <= t0_i \
                            and not np.allclose(SSI_root_array[0], t0_i)) \
                        or (SSI_root_array[0] > t1_i \
                            and not np.allclose(SSI_root_array[0], t1_i)):
                        # Suppose we are in the ]t1_i,TT_col_nff_i] phase.
                        GL = GL_3
                        SSI = HS - GL
                        # Find (SSI - j) real root.
                        SSI_root_array = get_real_roots(SSI - j)
                        if SSI_root_array.size == 0 \
                            or (SSI_root_array[0] <= t1_i \
                                and not np.allclose(SSI_root_array[0], t1_i)) \
                            or (SSI_root_array[0] > TT_col_nff_i \
                                and not np.allclose(SSI_root_array[0], TT_col_nff_i)):
                            # We must be in the ]TT_col_nff_i,infinity[ phase.
                            GL = GL_4
                            SSI = HS - GL
                            # Find (SSI - j) real root.
                            SSI_root_array = get_real_roots(SSI - j)
                            if SSI_root_array.size == 0 \
                                or (SSI_root_array[0] <= TT_col_nff_i \
                                    and not np.allclose(SSI_root_array[0], TT_col_nff_i)):
                                raise Exception('ERROR !!!!! This shouldn\'t occurred')
                            if HS(SSI_root_array[0]) > Nff_i:
                                HS = np.poly1d([Nff_i])
                            if GL(SSI_root_array[0]) < 0.0:
                                GL = np.poly1d([0.0])
                            SSI = HS - GL
                            # Find (SSI - j) real root again.
                            SSI_root_array = get_real_roots(SSI - j)
                            if SSI_root_array.size == 0 \
                                or (SSI_root_array[0] <= TT_col_nff_i \
                                    and not np.allclose(SSI_root_array[0], TT_col_nff_i)):
                                raise Exception('ERROR !!!!! This shouldn\'t occurred')   
                TT_sen_phytomer_i_list.append(SSI_root_array[0])
        
        TT_sen_phytomer_list += TT_sen_phytomer_i_list
    
    return TT_sen_phytomer_list


def _gen_absolute_TT_del_phytomer_list(id_phen_list, TT_sen_phytomer_list, dynT_dataframe, decimal_elongated_internode_number):
    '''Generate the *TT_del_phytomer* column.'''
    TT_del_phytomer_series = pandas.Series(np.nan, index=range(len(id_phen_list)))
    df = pandas.DataFrame({'id_phen': id_phen_list, 'TT_sen_phytomer': TT_sen_phytomer_list})
    Nbr_Fhaut_persistant = round(decimal_elongated_internode_number) + 1
    for name, group in df.groupby('id_phen'):
        dynT_dataframe_i = dynT_dataframe[dynT_dataframe['id_axis'] == name]
        a_cohort_i = dynT_dataframe_i['a_cohort'][dynT_dataframe_i.first_valid_index()]
        Nff_i = dynT_dataframe_i['Nff'][dynT_dataframe_i.first_valid_index()]
        N_cohort_i = dynT_dataframe_i['N_cohort'][dynT_dataframe_i.first_valid_index()]
        TT_del_Fhaut_phytomer_index = Nff_i - Nbr_Fhaut_persistant + N_cohort_i
        TT_del_phytomer_series[group.index[:TT_del_Fhaut_phytomer_index]] = group[:TT_del_Fhaut_phytomer_index]['TT_sen_phytomer'] + params.DELAIS_PHYLL_SEN_DISP / a_cohort_i
        TT_del_phytomer_series[group.index[TT_del_Fhaut_phytomer_index:]] = params.TT_DEL_FHAUT
    return TT_del_phytomer_series.tolist()
    
    
def create_phenT_first(phenT_abs_dataframe):
    '''
    Create the :ref:`phenT_first <phenT_first>` dataframe.

    :Parameters:
    
        - `phenT_abs_dataframe` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
        
    :Returns:
        The :ref:`phenT_first <phenT_first>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * *phenT_abs_dataframe* must be completely filled, i.e. must not contain 
          any NA value.
    
    '''
    tools.checkValidity(phenT_abs_dataframe.count().max() == phenT_abs_dataframe.count().min() == phenT_abs_dataframe.index.size)
    # Create a dataframe for first leaf (i.e. index_phytomer == 1) from phenT_abs_dataframe
    def first_leaf_criterion(index_i):
        # Permits to select first leaves row (i.e. row for which index_phytomer == 1).
        if phenT_abs_dataframe['index_phytomer'][index_i] == 1:
            return True
        return False
    phenT_first_dataframe = phenT_abs_dataframe.select(first_leaf_criterion)
    return phenT_first_dataframe


def create_phenT(phenT_abs_dataframe, phenT_first_dataframe):
    '''
    Create the :ref:`phenT <phenT>` dataframe.

    :Parameters:
    
        - `phenT_abs_dataframe` (:class:`pandas.DataFrame`) - the :ref:`phenT_abs <phenT_abs>` dataframe.
        - `phenT_first_dataframe` (:class:`pandas.DataFrame`) - the :ref:`phenT_first <phenT_first>` dataframe.
        
    :Returns:
        The :ref:`phenT <phenT>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    .. warning:: 
    
        * *phenT_abs_dataframe* must be completely filled, i.e. must not contain 
          any NA value.
    
    '''
    tools.checkValidity(phenT_abs_dataframe.count().max() == phenT_abs_dataframe.count().min() == phenT_abs_dataframe.index.size)
    phenT_dataframe = pandas.DataFrame(index=phenT_abs_dataframe.index, columns=['id_phen', 'index_rel_phytomer', 'dTT_em_phytomer', 'dTT_col_phytomer', 'dTT_sen_phytomer', 'dTT_del_phytomer'], dtype=float)
    phenT_dataframe['id_phen'] = phenT_abs_dataframe['id_phen']
    tmp_series = pandas.Series(phenT_dataframe.index)
    for name, group in phenT_abs_dataframe.groupby('id_phen'):
        current_first_leaf_row = phenT_first_dataframe[phenT_first_dataframe['id_phen'] == name]
        def normalize_index_phytomer(i):
            return group['index_phytomer'][i] / float(group['index_phytomer'][group.index[-1]])
        phenT_dataframe['index_rel_phytomer'][group.index] = tmp_series[group.index].map(normalize_index_phytomer)
        def get_relative_TT_col_phytomer(i):
            return group['TT_col_phytomer'][i] - current_first_leaf_row['TT_col_phytomer'][current_first_leaf_row.first_valid_index()]
        phenT_dataframe['dTT_col_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_col_phytomer)
        def get_relative_TT_em_phytomer(i):
            return group['TT_em_phytomer'][i] - current_first_leaf_row['TT_em_phytomer'][current_first_leaf_row.first_valid_index()]
        phenT_dataframe['dTT_em_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_em_phytomer)
        def get_relative_TT_sen_phytomer(i):
            return group['TT_sen_phytomer'][i] - current_first_leaf_row['TT_sen_phytomer'][current_first_leaf_row.first_valid_index()]
        phenT_dataframe['dTT_sen_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_sen_phytomer)
        def get_relative_TT_del_phytomer(i):
            return group['TT_del_phytomer'][i] - current_first_leaf_row['TT_del_phytomer'][current_first_leaf_row.first_valid_index()]
        phenT_dataframe['dTT_del_phytomer'][group.index] = tmp_series[group.index].map(get_relative_TT_del_phytomer)
        
    return phenT_dataframe


def _gen_HS_GL_polynomial(HS_break_i, N_cohort_i, id_axis_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, n0_i, n1_i, n2_i, t0_i, t1_i, a_i, c_i):
    MS = N_cohort_i == 1.0
    a2_i = (Nff_i - HS_break_i) / (TT_col_nff_i - TT_col_break_i)
    # define HS(TT)
    HS_1 = np.poly1d([a_cohort_i, - a_cohort_i * TT_col_0_i]) # j < HS_break_i
    HS_2 = np.poly1d([a2_i, - a2_i * TT_col_break_i + HS_break_i]) # j >= HS_break_i
    # define GL(TT) for all phases except TT < t0_i (because it depends on j)
    if MS:
        GL_2 = np.poly1d([(n1_i - n0_i) / (t1_i - t0_i), n0_i - t0_i * (n1_i - n0_i) / (t1_i - t0_i)])
        GL_3 = np.poly1d([(n2_i - n1_i) / (TT_col_nff_i - t1_i), n1_i - t1_i * (n2_i - n1_i) / (TT_col_nff_i - t1_i)])
    else: # tillers
        GL_2 = np.poly1d([n0_i])
        GL_3 = np.poly1d([(n2_i - n0_i) / (TT_col_nff_i - t1_i), n0_i - t1_i * (n2_i - n0_i) / (TT_col_nff_i - t1_i)])
    GL_4 = np.poly1d([a_i, - 3 * a_i * TT_col_nff_i, 3 * a_i * TT_col_nff_i**2 + c_i, - a_i * TT_col_nff_i**3 - c_i * TT_col_nff_i + n2_i])
    return HS_1, HS_2, GL_2, GL_3, GL_4


def create_HS_GL_SSI_T(dynT_dataframe):
    '''
    Create the :ref:`HS_GL_SSI_T <HS_GL_SSI_T>` dataframe.

    :Parameters:
    
        - `dynT_dataframe` (:class:`pandas.DataFrame`) - the :ref:`dynT <dynT>` dataframe.
        
    :Returns:
        The :ref:`HS_GL_SSI_T <HS_GL_SSI_T>` dataframe.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    '''
    HS_GL_SSI_dynamic_dataframe = pandas.DataFrame(columns=['id_axis', 'TT', 'HS', 'GL', 'SSI'])
    for i in dynT_dataframe.index:
        N_cohort_i, id_axis_i, cardinality_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, dTT_MS_cohort_i, n0_i, n1_i, n2_i, t0_i, t1_i, hs_t1_i, a_i, c_i, RMSE_gl = dynT_dataframe.ix[i].tolist()
        HS_break_i = a_cohort_i * (TT_col_break_i - TT_col_0_i)
        HS_1, HS_2, GL_2, GL_3, GL_4 = _gen_HS_GL_polynomial(HS_break_i, N_cohort_i, id_axis_i, Nff_i, a_cohort_i, TT_col_0_i, TT_col_break_i, TT_col_nff_i, n0_i, n1_i, n2_i, t0_i, t1_i, a_i, c_i)
        
        t0_i, t1_i, TT_col_nff_i, TT_col_break_i = np.round([t0_i, t1_i, TT_col_nff_i, TT_col_break_i]).astype(int)
        
        TT_1_1 = np.arange(0, t0_i)
        TT_1_2 = np.arange(t0_i, t0_i)
        TT_2_1 = np.arange(t0_i, t1_i)
        TT_2_2 = np.arange(t1_i, t1_i)
        TT_3_1 = np.arange(t1_i, TT_col_nff_i)
        TT_3_2 = np.arange(TT_col_nff_i, TT_col_nff_i)
        TT_4_1 = np.arange(TT_col_nff_i, params.TT_DEL_FHAUT)
        TT_4_2 = np.arange(params.TT_DEL_FHAUT, params.TT_DEL_FHAUT)
        
        if TT_col_break_i != 0.0: # bilinear mode
            if TT_col_break_i <= t0_i:
                TT_1_1 = np.arange(0, TT_col_break_i)
                TT_1_2 = np.arange(TT_col_break_i, t0_i)
            elif TT_col_break_i <= t1_i:
                TT_2_1 = np.arange(t0_i, TT_col_break_i)
                TT_2_2 = np.arange(TT_col_break_i, t1_i)
            elif TT_col_break_i <= TT_col_nff_i:
                TT_3_1 = np.arange(t1_i, TT_col_break_i)
                TT_3_2 = np.arange(TT_col_break_i, TT_col_nff_i)
            else:
                TT_4_1 = np.arange(TT_col_nff_i, TT_col_break_i)
                TT_4_2 = np.arange(TT_col_break_i, params.TT_DEL_FHAUT)
            
        HS_1_TT_1_1 = np.clip(HS_1(TT_1_1), 0.0, Nff_i)
        HS_1_TT_2_1 = np.clip(HS_1(TT_2_1), 0.0, Nff_i)
        HS_1_TT_3_1 = np.clip(HS_1(TT_3_1), 0.0, Nff_i)
        HS_1_TT_4_1 = np.clip(HS_1(TT_4_1), 0.0, Nff_i)
        HS_1_TT_1_2 = np.clip(HS_1(TT_1_2), 0.0, Nff_i)
        HS_1_TT_2_2 = np.clip(HS_1(TT_2_2), 0.0, Nff_i)
        HS_1_TT_3_2 = np.clip(HS_1(TT_3_2), 0.0, Nff_i)
        HS_1_TT_4_2 = np.clip(HS_1(TT_4_2), 0.0, Nff_i)
        
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
        
        HS_GL_SSI_dynamic_dataframe_i = pandas.DataFrame(index=np.arange(params.TT_DEL_FHAUT), columns=HS_GL_SSI_dynamic_dataframe.columns)
        HS_GL_SSI_dynamic_dataframe_i['TT'] = pandas.Series(np.concatenate((TT_1_1, TT_1_2, TT_2_1, TT_2_2, TT_3_1, TT_3_2, TT_4_1, TT_4_2)))
        HS_GL_SSI_dynamic_dataframe_i['HS'] = pandas.Series(np.concatenate((HS_1_TT_1_1, HS_1_TT_1_2, HS_1_TT_2_1, HS_1_TT_2_2, HS_1_TT_3_1, HS_1_TT_3_2, HS_1_TT_4_1, HS_1_TT_4_2)))
        HS_GL_SSI_dynamic_dataframe_i['GL'] = pandas.Series(np.concatenate((GL_1_TT_1_1, GL_1_TT_1_2, GL_2_TT_2_1, GL_2_TT_2_2, GL_3_TT_3_1, GL_3_TT_3_2, GL_4_TT_4_1, GL_4_TT_4_2)))
        HS_GL_SSI_dynamic_dataframe_i['SSI'] = pandas.Series(np.concatenate((SSI_1_TT_1_1, SSI_1_TT_1_2, SSI_2_TT_2_1, SSI_2_TT_2_2, SSI_3_TT_3_1, SSI_3_TT_3_2, SSI_4_TT_4_1, SSI_4_TT_4_2)))
        HS_GL_SSI_dynamic_dataframe_i['id_axis'] = np.repeat(id_axis_i, HS_GL_SSI_dynamic_dataframe_i.index.size)
        HS_GL_SSI_dynamic_dataframe = HS_GL_SSI_dynamic_dataframe.append(HS_GL_SSI_dynamic_dataframe_i, ignore_index=True)
    
    return HS_GL_SSI_dynamic_dataframe


    