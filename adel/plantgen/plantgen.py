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
Front end for the generation of the input data expected by ADEL. User should 
look at this module first. One can then look at the other modules of :mod:`alinea.adel.plantgen` 
for additional information. 

Authors: Mariem Abichou, Camille Chambon, Bruno Andrieu 
'''

import numpy as np
import pandas
import warnings

from adel.plantgen import axeT, dimT, dynT, phenT, tools, params

warnings.simplefilter('always', tools.InputWarning)

class DataCompleteness:
    '''
    This enumerate defines the different degrees of completeness that the data 
    documented by the user can have.
    
    .. seealso:: :func:`gen_adel_input_data`  
    '''
    MIN='MIN'
    SHORT='SHORT'
    FULL='FULL'


def gen_adel_input_data(dynT_user,
                        dimT_user,
                        plants_number=100, 
                        plants_density=250,
                        decide_child_axis_probabilities={'T0': 0.0, 'T1': 0.900, 'T2': 0.983, 'T3': 0.817, 'T4': 0.117}, 
                        MS_leaves_number_probabilities={'10': 0.145, '11': 0.818, '12': 0.037, '13': 0.0, '14': 0.0},
                        ears_density=500,
                        GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                        delais_TT_stop_del_axis=600,
                        TT_col_break=0.0,
                        inner_params={}                        
                        ):
    '''
    Create the dataframes which contain the plant data to be used as input for 
    generating plot with ADEL, and some other dataframes for debugging purpose.
    Also create a dictionary which stores the values of the arguments of 
    :func:`gen_adel_input_data <alinea.adel.plantgen.plantgen.gen_adel_input_data>`. 
    This dictionary is aimed to log the configuration used for the construction.
    
    See :ref:`adel_input` for a description of the input tables expected by ADEL, 
    and :ref:`plantgen` for a description of the dataframes created for debug. 
    
    Different degrees of completeness of data provided by the user are acceptable. 
    See :ref:`levels_of_completeness` for more details.
    
    The dataframes are created as follows:
        * initialization of the following dataframes:
            * *axeT_tmp*, calling :func:`alinea.adel.plantgen.axeT.create_axeT_tmp`,
            * *dynT_tmp*, calling :func:`alinea.adel.plantgen.dynT.create_dynT_tmp`.
            * *dimT_tmp*, calling :func:`alinea.adel.plantgen.dimT.create_dimT_tmp`
            * :ref:`tilleringT <tilleringT>`, calling :func:`alinea.adel.plantgen.axeT.create_tilleringT`
            * :ref:`cardinalityT <cardinalityT>`, calling :func:`alinea.adel.plantgen.axeT.create_cardinalityT`
        * filling of the dataframes *dynT_user* and *dimT_user* set by the user 
        * calculate the number of elongated internodes
        * construction of the following dataframes:
            * :ref:`dynT <dynT>`, calling :func:`alinea.adel.plantgen.dynT.create_dynT`,
            * :ref:`phenT_abs <phenT_abs>`, calling :func:`alinea.adel.plantgen.phenT.create_phenT_abs`,
            * :ref:`phenT_first <phenT_first>`, calling :func:`alinea.adel.plantgen.phenT.create_phenT_first`,
            * :ref:`phenT <phenT>`, calling :func:`alinea.adel.plantgen.phenT.create_phenT`,
            * :ref:`axeT <axeT>`, calling :func:`alinea.adel.plantgen.axeT.create_axeT`,
            * :ref:`dimT_abs <dimT_abs>`, calling :func:`alinea.adel.plantgen.dimT.create_dimT_abs`,
            * :ref:`dimT <dimT>`, calling :func:`alinea.adel.plantgen.dimT.create_dimT`,
            * :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`, calling :func:`alinea.adel.plantgen.phenT.create_HS_GL_SSI_T`,

        These tables are returned to be used as ADEL input:
            * the :ref:`axeT <axeT>`, 
            * the :ref:`dimT <dimT>`, 
            * the :ref:`phenT <phenT>`.
          
        These tables are intermediate tables and are returned for debugging purpose:
            * the :ref:`tilleringT <tilleringT>`,
            * the :ref:`cardinalityT <cardinalityT>`,
            * the :ref:`phenT_abs <phenT_abs>`,
            * the :ref:`dimT_abs <dimT_abs>`,
            * the :ref:`dynT <dynT>`, 
            * the :ref:`phenT_first <phenT_first>`,
            * the :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`,
        
    :Parameters:
    
        - `dynT_user` (:class:`pandas.DataFrame`) - the leaf dynamic 
          parameters set by the user.
                
        - `dimT_user` (:class:`pandas.DataFrame`) - the dimensions of the organs 
          set by the user. 
              
        - `plants_number` (:class:`int`) - the number of plants to be generated.
        
        - `plants_density` (:class:`int`) - the number of plants that are present 
          after loss due to bad emergence, early death..., per square meter.
        
        - `decide_child_axis_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          for each child cohort the probability of emergence of an axis when the parent 
          axis is present. The keys are the identifiers of the child axes ('T0', 'T1', 
          'T2', ...) and the values are the probabilities.
        
        - `MS_leaves_number_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          the probability distribution of the final number of main stem leaves. 
          The keys are the final numbers of main stem leaves, and the values are 
          the probabilities distribution.
          
        - `ears_density` (:class:`int`) - the number of ears per square meter. 
        
        - `GL_number` (:class:`dict` of :class:`float`::class:`float`) - the GL decimal numbers measured at 
          several thermal times (including the senescence end). The keys are the 
          thermal times, and the values are the GL decimal numbers.
        
        - `delais_TT_stop_del_axis` (:class:`int`) - This variable represents the time in 
          thermal time between an axis stop growing and its disappearance (it 
          concerns only the axes that do not regress and which do not produce any 
          cob).
        
        - `TT_col_break` (:class:`float`) - the thermal time when the rate of Haun Stage 
          is changing. If phyllochron is constant, then *TT_col_break* is null.
          
        - `inner_params` (:class:`dict`) - the values of the inner parameters used 
          for the construction of the input tables. These parameters are the same 
          as the ones defined in the module :mod:`params <alinea.adel.plantgen.params>`. 
          *inner_params* is NOT mandatory: if not all inner parameters are documented 
          in *inner_params*, then we use the default values defined in :mod:`params <alinea.adel.plantgen.params>` 
          for the inner parameters which are missing.   
        
    :Returns:
        Return :ref:`axeT <axeT>`, :ref:`dimT <dimT>`, 
        :ref:`phenT <phenT>`, :ref:`phenT_abs <phenT_abs>`, :ref:`dimT_abs <dimT_abs>`, 
        :ref:`dynT <dynT>`, :ref:`phenT_first <phenT_first>`, :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`, 
        :ref:`tilleringT <tilleringT>`, :ref:`cardinalityT <cardinalityT>`, and 
        a dictionary which stores the configuration used for the construction.
    
    :Returns Type:
        tuple
        
    .. seealso:: :class:`DataCompleteness`
                 :mod:`alinea.adel.plantgen.axeT`
                 :mod:`alinea.adel.plantgen.dimT`
                 :mod:`alinea.adel.plantgen.dynT`
                 :mod:`alinea.adel.plantgen.params`
                 :mod:`alinea.adel.plantgen.phenT`
                 :mod:`alinea.adel.plantgen.tools`
                 
    '''
    # save the name and the value of each argument
    config = locals()
    
    if sum(MS_leaves_number_probabilities.values()) != 1.0:
        raise tools.InputError("the sum of the probabilities defined in MS_leaves_number_probabilities is not equal to 1.0")
    
    # update values defined in alinea.adel.plantgen.params from values in inner_params
    attribute_names = set(dir(params))
    attribute_names.intersection_update(inner_params.keys())
    params.__dict__.update(dict([(key, value) for key, value in inner_params.items() if key in attribute_names]))
    
    possible_axes = \
        set([id_axis for (id_axis, probability) in
             decide_child_axis_probabilities.iteritems() if probability != 0.0])
        
    possible_MS_N_phytomer_potential = \
        set([MS_N_phytomer_potential for (MS_N_phytomer_potential, probability) in
             MS_leaves_number_probabilities.iteritems() if probability != 0.0])
    
    # check plants_number, decide_child_axis_probabilities, plants_density and ears_density validity
    decide_child_cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(decide_child_axis_probabilities)
    (theoretical_cohort_cardinalities, 
     theoretical_axis_cardinalities) = tools.calculate_theoretical_cardinalities(plants_number, 
                                                                                 decide_child_cohort_probabilities,
                                                                                 decide_child_axis_probabilities,
                                                                                 params.FIRST_CHILD_DELAY)
    theoretical_cardinalities_sum = sum(theoretical_cohort_cardinalities.values())
    number_of_ears = plants_number * ears_density / float(plants_density)
    if number_of_ears < plants_number:
        raise tools.InputError("the number of ears (%s) is lesser than plants_number (%s). \
The number of ears (%s) is calculated from plants_number (%s), ears_density (%s) \
and plants_density (%s)." % (number_of_ears, plants_number, number_of_ears, plants_number, 
                             ears_density, plants_density))
    
    if number_of_ears > theoretical_cardinalities_sum:
        raise tools.InputError("the number of ears (%s) is greater than the theoretical number of axes (%s). \
The number of ears (%s) is calculated from plants_number (%s), ears_density (%s) and plants_density (%s). \
The theoretical number of axes (%s) is calculated from plants_number (%s), decide_child_cohort_probabilities \
(%s) and params.FIRST_CHILD_DELAY (%s)" % (number_of_ears, theoretical_cardinalities_sum, number_of_ears, 
                                           plants_number, ears_density, plants_density, theoretical_cardinalities_sum,
                                           plants_number, decide_child_cohort_probabilities, params.FIRST_CHILD_DELAY))
    
    available_axes_warning_message = "the probabilities defined in decide_child_axis_probabilities (%s) and \
the axes documented by the user (%s) in %s indicate that some of the possible axes (%s) \
are not documented by the user. After the generation of the axes, if not all generated axes are documented by the user, \
then this will lead to an error."

    available_MS_N_phytomer_potential_warning_message = "the probabilities defined in MS_leaves_number_probabilities (%s) and \
the N_phytomer_potential of the MS documented by the user (%s) in %s indicate that some of the N_phytomer_potential of the MS (%s) \
are not documented by the user. After the generation of the phytomers of the MS, if not all generated phytomers \
of the MS are documented by the user, then this will lead to an error."
    
    # calculate dynT_user completeness and check its validity
    if 'N_phytomer_potential' in dynT_user.columns:
        dynT_user_completeness = DataCompleteness.FULL
        expected_dynT_user_columns = ['id_axis', 'N_phytomer_potential', 'a_cohort', 'TT_col_0', 'TT_col_N_phytomer_potential', 'n0', 'n1', 'n2']
        if dynT_user.columns.tolist() != expected_dynT_user_columns:
            raise tools.InputError("dynT_user does not have the columns: %s" % ', '.join(expected_dynT_user_columns))
        grouped = dynT_user.groupby(['id_axis', 'N_phytomer_potential'])
        if len(grouped.groups) != dynT_user.index.size:
            raise tools.InputError("dynT_user contains duplicated (id_axis, N_phytomer_potential) pair(s)")
        available_axes = set(dynT_user['id_axis'].tolist())
        if not possible_axes.issubset(available_axes):
            warnings.warn(available_axes_warning_message % (decide_child_axis_probabilities,
                                                            list(available_axes),
                                                            'dynT_user',
                                                            list(possible_axes)),
                          tools.InputWarning)
        available_MS_N_phytomer_potential = set(dynT_user[dynT_user['id_axis'] == 'MS']['N_phytomer_potential'].tolist())
        if not possible_MS_N_phytomer_potential.issubset(available_MS_N_phytomer_potential):
            warnings.warn(available_MS_N_phytomer_potential_warning_message % (MS_leaves_number_probabilities,
                                                                     list(available_MS_N_phytomer_potential),
                                                                     'dynT_user',
                                                                     list(possible_MS_N_phytomer_potential)),
                          tools.InputWarning)
            
    elif dynT_user.count().max() == dynT_user.count().min() == dynT_user.index.size:
        dynT_user_completeness = DataCompleteness.SHORT
        expected_dynT_user_columns = ['id_axis', 'a_cohort', 'TT_col_0', 'TT_col_N_phytomer_potential', 'n0', 'n1', 'n2']
        if dynT_user.columns.tolist() != expected_dynT_user_columns:
            raise tools.InputError("dynT_user does not have the columns: %s" % ', '.join(expected_dynT_user_columns))
        if dynT_user['id_axis'].unique().size != dynT_user['id_axis'].size:
            raise tools.InputError("dynT_user contains duplicated id_axis")
        available_axes = set(dynT_user['id_axis'].tolist())
        if not possible_axes.issubset(available_axes):
            warnings.warn(available_axes_warning_message % (decide_child_axis_probabilities,
                                                            list(available_axes),
                                                            'dynT_user',
                                                            list(possible_axes)),
                          tools.InputWarning)
        
    else:
        dynT_user_completeness = DataCompleteness.MIN
        expected_dynT_user_columns = ['id_axis', 'a_cohort', 'TT_col_0', 'TT_col_N_phytomer_potential', 'n0', 'n1', 'n2']
        if dynT_user.columns.tolist() != expected_dynT_user_columns:
            raise tools.InputError("dynT_user does not have the columns: %s" % ', '.join(expected_dynT_user_columns))
        available_axes = set(dynT_user['id_axis'].tolist())
        if not possible_axes.issubset(available_axes):
            warnings.warn(available_axes_warning_message % (decide_child_axis_probabilities,
                                                            list(available_axes),
                                                            'dynT_user',
                                                            list(possible_axes)), 
                          tools.InputWarning)
        
    # calculate dimT_user completeness and check its validity
    if 'N_phytomer_potential' in dimT_user.columns:
        dimT_user_completeness = DataCompleteness.FULL
        expected_dimT_user_columns = ['id_axis', 'N_phytomer_potential', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        if dimT_user.columns.tolist() != expected_dimT_user_columns:
            raise tools.InputError("dimT_user does not have the columns: %s" % ', '.join(expected_dimT_user_columns))
        grouped = dimT_user.groupby(['id_axis', 'N_phytomer_potential', 'index_phytomer'])
        if len(grouped.groups) != dimT_user.index.size:
            raise tools.InputError("dimT_user contains duplicated (id_axis, N_phytomer_potential, index_phytomer) triplet(s)")
        available_axes = set(dimT_user['id_axis'].tolist())
        if not possible_axes.issubset(available_axes):
            warnings.warn(available_axes_warning_message % (decide_child_axis_probabilities,
                                                            list(available_axes),
                                                            'dimT_user',
                                                            list(possible_axes)),
                          tools.InputWarning)
        available_MS_N_phytomer_potential = set(dimT_user[dimT_user['id_axis'] == 'MS']['N_phytomer_potential'].tolist())
        if not possible_MS_N_phytomer_potential.issubset(available_MS_N_phytomer_potential):
            warnings.warn(available_MS_N_phytomer_potential_warning_message % (MS_leaves_number_probabilities,
                                                                     list(available_MS_N_phytomer_potential),
                                                                     'dimT_user',
                                                                     list(possible_MS_N_phytomer_potential)),
                          tools.InputWarning)
        
    elif 'id_axis' in dimT_user.columns:
        dimT_user_completeness = DataCompleteness.SHORT
        expected_dimT_user_columns = ['id_axis', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        if dimT_user.columns.tolist() != expected_dimT_user_columns:
            raise tools.InputError("dimT_user does not have the columns: %s" % ', '.join(expected_dimT_user_columns))
        grouped = dimT_user.groupby(['id_axis', 'index_phytomer'])
        if len(grouped.groups) != dimT_user.index.size:
            raise tools.InputError("dimT_user contains duplicated (id_axis, index_phytomer) pair(s)")
        available_axes = set(dimT_user['id_axis'].tolist())
        if not possible_axes.issubset(available_axes):
            warnings.warn(available_axes_warning_message % (decide_child_axis_probabilities,
                                                            list(available_axes),
                                                            'dimT_user',
                                                            list(possible_axes)),
                          tools.InputWarning)
        max_available_MS_N_phytomer_potential = dimT_user[dimT_user['id_axis'] == 'MS']['index_phytomer'].max()
        if max(possible_MS_N_phytomer_potential) > max_available_MS_N_phytomer_potential:
            warnings.warn(available_MS_N_phytomer_potential_warning_message % (MS_leaves_number_probabilities,
                                                                     ', '.join([str(max_available_MS_N_phytomer_potential)]),
                                                                     'dimT_user',
                                                                     ', '.join([str(max(possible_MS_N_phytomer_potential))])),
                          tools.InputWarning)
        
    else:
        dimT_user_completeness = DataCompleteness.MIN
        expected_dimT_user_columns = ['index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        if dimT_user.columns.tolist() != expected_dimT_user_columns:
            raise tools.InputError("dimT_user does not have the columns: %s" % ', '.join(expected_dimT_user_columns))
        if dimT_user['index_phytomer'].unique().size != dimT_user['index_phytomer'].size:
            raise tools.InputError("dimT_user contains duplicated index_phytomer")
        max_available_MS_N_phytomer_potential = dimT_user['index_phytomer'].max()
        if max(possible_MS_N_phytomer_potential) > max_available_MS_N_phytomer_potential:
            warnings.warn(available_MS_N_phytomer_potential_warning_message % (MS_leaves_number_probabilities,
                                                                     ', '.join([str(max_available_MS_N_phytomer_potential)]),
                                                                     'dimT_user',
                                                                     ', '.join([str(max(possible_MS_N_phytomer_potential))])),
                          tools.InputWarning)
            
    
    # 2. first step of the fit process
    (axeT_tmp, 
     dimT_tmp, 
     dynT_tmp, 
     cardinalityT) = _gen_adel_input_data_first(plants_number, 
                                                decide_child_cohort_probabilities, 
                                                MS_leaves_number_probabilities, 
                                                theoretical_cohort_cardinalities,
                                                theoretical_axis_cardinalities)
     
    most_frequent_dynT_tmp = pandas.DataFrame(columns=dynT_tmp.columns)
    for id_axis, dynT_tmp_group in dynT_tmp.groupby('id_axis'):
        idxmax = dynT_tmp_group['cardinality'].idxmax()
        most_frequent_dynT_tmp = pandas.concat([most_frequent_dynT_tmp, dynT_tmp_group.ix[idxmax:idxmax]], ignore_index=True)
    most_frequent_dynT_tmp_grouped = most_frequent_dynT_tmp.groupby(['id_axis', 'N_phytomer_potential'])
    
    # 3. complete dynT_tmp
    if dynT_user_completeness == DataCompleteness.MIN:
        dynT_user_grouped = dynT_user.groupby('id_axis')
        MS_dynT_user = dynT_user_grouped.get_group('MS')
        MS_TT_col_N_phytomer_potential = MS_dynT_user['TT_col_N_phytomer_potential'][MS_dynT_user.first_valid_index()]
        for id_axis, dynT_tmp_group in dynT_tmp.groupby('id_axis'):
            if not dynT_user['id_axis'].isin([id_axis]).any():
                id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
            dynT_user_group = dynT_user_grouped.get_group(id_axis)
            index_to_get = dynT_user_group.index[0]
            index_to_set = dynT_tmp_group.index[0]
            if id_axis == 'MS':
                for column in ['a_cohort', 'TT_col_0', 'TT_col_N_phytomer_potential', 'n0', 'n1', 'n2']:
                    dynT_tmp[column][index_to_set] = dynT_user[column][index_to_get]
            dynT_tmp['TT_col_break'][index_to_set] = TT_col_break
            current_TT_col_N_phytomer_potential = dynT_user_group['TT_col_N_phytomer_potential'][index_to_get]
            current_dTT_MS_cohort = current_TT_col_N_phytomer_potential - MS_TT_col_N_phytomer_potential
            dynT_tmp['dTT_MS_cohort'][index_to_set] = current_dTT_MS_cohort
    elif dynT_user_completeness == DataCompleteness.SHORT:
        dynT_user_grouped = dynT_user.groupby('id_axis')
        MS_dynT_user = dynT_user_grouped.get_group('MS')
        MS_TT_col_N_phytomer_potential = MS_dynT_user['TT_col_N_phytomer_potential'][MS_dynT_user.first_valid_index()]
        for id_axis, dynT_tmp_group in dynT_tmp.groupby('id_axis'):
            if not dynT_user['id_axis'].isin([id_axis]).any():
                id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
            dynT_user_group = dynT_user_grouped.get_group(id_axis)
            index_to_get = dynT_user_group.index[0]
            index_to_set = dynT_tmp_group.index[0]
            for column in ['a_cohort', 'TT_col_0', 'TT_col_N_phytomer_potential', 'n0', 'n1', 'n2']:
                dynT_tmp[column][index_to_set] = dynT_user[column][index_to_get]
            dynT_tmp['TT_col_break'][index_to_set] = TT_col_break
            current_TT_col_N_phytomer_potential = dynT_user_group['TT_col_N_phytomer_potential'][index_to_get]
            current_dTT_MS_cohort = current_TT_col_N_phytomer_potential - MS_TT_col_N_phytomer_potential
            dynT_tmp['dTT_MS_cohort'][index_to_set] = current_dTT_MS_cohort
    elif dynT_user_completeness == DataCompleteness.FULL:
        dynT_user_grouped = dynT_user.groupby(['id_axis', 'N_phytomer_potential'])
        MS_N_phytomer_potential = dynT_tmp['N_phytomer_potential'][dynT_tmp[dynT_tmp['id_axis'] == 'MS']['cardinality'].idxmax()]
        MS_dynT_user = dynT_user_grouped.get_group(('MS', MS_N_phytomer_potential))
        MS_TT_col_N_phytomer_potential = MS_dynT_user['TT_col_N_phytomer_potential'][MS_dynT_user.first_valid_index()]
        for (id_axis, N_phytomer_potential), dynT_tmp_group in dynT_tmp.groupby(['id_axis', 'N_phytomer_potential']):
            if not dynT_user['id_axis'].isin([id_axis]).any():
                if most_frequent_dynT_tmp_grouped.groups.has_key((id_axis, N_phytomer_potential)):
                    id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
                else:
                    continue
            if not dynT_user_grouped.groups.has_key((id_axis, N_phytomer_potential)):
                if most_frequent_dynT_tmp_grouped.groups.has_key((id_axis, N_phytomer_potential)):
                    raise tools.InputError("Dynamic of %s not documented" % ((id_axis, N_phytomer_potential),))
                else:
                    most_frequent_dynT_tmp_id_axis = most_frequent_dynT_tmp[most_frequent_dynT_tmp['id_axis'] == id_axis]
                    N_phytomer_potential = most_frequent_dynT_tmp_id_axis['N_phytomer_potential'][most_frequent_dynT_tmp_id_axis.first_valid_index()]
                    if not dynT_user_grouped.groups.has_key((id_axis, N_phytomer_potential)):
                        raise tools.InputError("Dynamic of %s not documented" % ((id_axis, N_phytomer_potential),))
            dynT_user_group = dynT_user_grouped.get_group((id_axis, N_phytomer_potential))
            index_to_get = dynT_user_group.index[0]
            index_to_set = dynT_tmp_group.index[0]
            for column in ['a_cohort', 'TT_col_0', 'TT_col_N_phytomer_potential', 'n0', 'n1', 'n2']:
                dynT_tmp[column][index_to_set] = dynT_user[column][index_to_get]
            dynT_tmp['TT_col_break'][index_to_set] = TT_col_break
            current_TT_col_N_phytomer_potential = dynT_user_group['TT_col_N_phytomer_potential'][index_to_get]
            current_dTT_MS_cohort = current_TT_col_N_phytomer_potential - MS_TT_col_N_phytomer_potential
            dynT_tmp['dTT_MS_cohort'][index_to_set] = current_dTT_MS_cohort

    # 4. complete dimT_tmp
    organ_dim_list = ['L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
    if dimT_user_completeness == DataCompleteness.MIN:
        MS_dynT_tmp = dynT_tmp[dynT_tmp['id_axis'] == 'MS']
        MS_most_frequent_N_phytomer_potential = MS_dynT_tmp['N_phytomer_potential'][MS_dynT_tmp['cardinality'].idxmax()]
        dimT_tmp_grouped = dimT_tmp.groupby(['id_axis', 'N_phytomer_potential'])
        dimT_tmp_indexes_to_set = dimT_tmp_grouped.groups[('MS', MS_most_frequent_N_phytomer_potential)]
        N_phytomer_potential_to_set = len(dimT_tmp_indexes_to_set)
        max_available_MS_N_phytomer_potential = dimT_user['index_phytomer'].max()
        if N_phytomer_potential_to_set > max_available_MS_N_phytomer_potential:
            raise tools.InputError("Dimensions of index_phytomer=%s not documented" % N_phytomer_potential_to_set)
        dimT_user_indexes_to_get = range(len(dimT_tmp_indexes_to_set))
        for organ_dim in organ_dim_list:
            dimT_tmp[organ_dim][dimT_tmp_indexes_to_set] = dimT_user[organ_dim][dimT_user_indexes_to_get]
    elif dimT_user_completeness == DataCompleteness.SHORT:
        dimT_user_grouped = dimT_user.groupby('id_axis')
        dynT_tmp_grouped = dynT_tmp.groupby('id_axis')
        for id_axis, dimT_tmp_group in dimT_tmp.groupby('id_axis'):
            dynT_tmp_group = dynT_tmp_grouped.get_group(id_axis)
            N_phytomer_potential = dynT_tmp_group['N_phytomer_potential'][dynT_tmp_group['cardinality'].idxmax()]
            indexes_to_set = dimT_tmp_group[dimT_tmp_group['N_phytomer_potential'] == N_phytomer_potential].index
            if not dimT_user['id_axis'].isin([id_axis]).any():
                id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
            dimT_user_group = dimT_user_grouped.get_group(id_axis)
            max_available_N_phytomer_potential = dimT_user_group['index_phytomer'].max()
            N_phytomer_potential_to_set = len(indexes_to_set)
            if N_phytomer_potential_to_set > max_available_N_phytomer_potential:
                raise tools.InputError("Dimensions of %s not documented" % ((id_axis, N_phytomer_potential_to_set),))
            indexes_to_get = dimT_user_group.index[:N_phytomer_potential_to_set]
            for organ_dim in organ_dim_list:
                dimT_tmp[organ_dim][indexes_to_set] = dimT_user[organ_dim][indexes_to_get]
    elif dimT_user_completeness == DataCompleteness.FULL:
        dimT_user_grouped = dimT_user.groupby(['id_axis', 'N_phytomer_potential'])
        for (id_axis, N_phytomer_potential), dimT_tmp_group in dimT_tmp.groupby(['id_axis', 'N_phytomer_potential']):
            if not dimT_user['id_axis'].isin([id_axis]).any():
                if most_frequent_dynT_tmp_grouped.groups.has_key((id_axis, N_phytomer_potential)):
                    id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
                else:
                    continue
            indexes_to_set = dimT_tmp_group.index
            if not dimT_user_grouped.groups.has_key((id_axis, N_phytomer_potential)):
                raise tools.InputError("Dimensions of %s not documented" % ((id_axis, N_phytomer_potential),))
            indexes_to_get = dimT_user_grouped.get_group((id_axis, N_phytomer_potential)).index
            for organ_dim in organ_dim_list:
                dimT_tmp[organ_dim][indexes_to_set] = dimT_user[organ_dim][indexes_to_get]
    
    # 5. second step of the fit process    
    (axeT_, 
     phenT_abs, 
     phenT_,
     dimT_abs, 
     dynT_, 
     phenT_first,
     HS_GL_SSI_T, 
     dimT_,
     tilleringT) = _gen_adel_input_data_second(axeT_tmp, 
                                               dimT_tmp, 
                                               dynT_tmp, 
                                               plants_number,
                                               plants_density,
                                               GL_number, 
                                               delais_TT_stop_del_axis, 
                                               ears_density,
                                               number_of_ears)
     
    return axeT_, dimT_, phenT_, phenT_abs, dimT_abs, dynT_, phenT_first, \
           HS_GL_SSI_T, tilleringT, cardinalityT, config


def _gen_adel_input_data_first(plants_number, 
                               decide_child_cohort_probabilities, 
                               MS_leaves_number_probabilities,
                               theoretical_cohort_cardinalities,
                               theoretical_axis_cardinalities):    
    '''Generate the input data: first step.'''
    # create axeT_tmp
    axeT_tmp = axeT.create_axeT_tmp(plants_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities)
    # create dynT_tmp
    dynT_tmp = dynT.create_dynT_tmp(axeT_tmp)
    # create dimT_tmp
    dimT_tmp = dimT.create_dimT_tmp(axeT_tmp)
    # create cardinalityT
    cardinalityT = axeT.create_cardinalityT(theoretical_cohort_cardinalities, 
                                            theoretical_axis_cardinalities, 
                                            axeT_tmp[['id_cohort', 'id_axis']])

    return axeT_tmp, dimT_tmp, dynT_tmp, cardinalityT


def _gen_adel_input_data_second(axeT_tmp, 
                                dimT_tmp, 
                                dynT_tmp,
                                plants_number, 
                                plants_density,
                                GL_number, 
                                delais_TT_stop_del_axis,
                                ears_density,
                                number_of_ears):
    '''Generate the input data: second step.'''
    # calculate decimal_elongated_internode_number
    decimal_elongated_internode_number = dynT.calculate_decimal_elongated_internode_number(dimT_tmp, dynT_tmp) 
    # create dynT
    dynT_ = dynT.create_dynT(dynT_tmp, GL_number, decimal_elongated_internode_number)
    # create phenT_tmp
    phenT_tmp = phenT.create_phenT_tmp(axeT_tmp, dynT_)
    # create phenT_first
    phenT_first = phenT.create_phenT_first(phenT_tmp)
    # calculate TT_flag_leaf_ligulation
    TT_flag_leaf_ligulation = dynT_['TT_col_N_phytomer_potential'][dynT_.first_valid_index()]
    # create axeT
    t1_most_frequent_MS = dynT_['t1'][dynT_.first_valid_index()]
    TT_regression_start = t1_most_frequent_MS + params.DELAIS_REG_MONT
    axeT_ = axeT.create_axeT(axeT_tmp, phenT_first, dynT_, TT_regression_start, TT_flag_leaf_ligulation, delais_TT_stop_del_axis, number_of_ears)
    # create tilleringT
    dynT_most_frequent_MS = dynT_.ix[dynT_.first_valid_index()]
    id_cohort_most_frequent_MS = str(dynT_most_frequent_MS['id_cohort'])
    N_phytomer_potential_most_frequent_MS = str(dynT_most_frequent_MS['N_phytomer_potential'])
    id_phen_most_frequent_MS = int(''.join([id_cohort_most_frequent_MS, N_phytomer_potential_most_frequent_MS]))
    TT_start = phenT_first['TT_app_phytomer'][phenT_first[phenT_first['id_phen'] == id_phen_most_frequent_MS].index[0]]
    tilleringT = axeT.create_tilleringT(TT_start, TT_regression_start, TT_flag_leaf_ligulation, plants_number, plants_density, axeT_tmp.index.size, ears_density)
    # create dimT_abs
    dimT_abs = dimT.create_dimT_abs(axeT_, dimT_tmp, phenT_tmp, dynT_)
    # create dimT
    dimT_ = dimT.create_dimT(dimT_abs)
    # create phenT_abs
    phenT_abs = phenT.create_phenT_abs(phenT_tmp, axeT_, dimT_abs)
    # create phenT
    phenT_ = phenT.create_phenT(phenT_abs, phenT_first)
    # create HS_GL_SSI_T 
    HS_GL_SSI_T = phenT.create_HS_GL_SSI_T(phenT_abs, axeT_tmp, dynT_)
    
    return axeT_, phenT_abs, phenT_, dimT_abs, dynT_, phenT_first, HS_GL_SSI_T, dimT_, tilleringT


def read_plantgen_inputs(inputs_filepath):
    '''
    Import the Python module at *inputs_filepath*, and return the args expected by 
    :func:`gen_adel_input_data`.
    
    :Parameters:
    
        - `inputs_filepath` (:class:`str`) - the Python module which contains the inputs 
          of :func:`gen_adel_input_data`.
          
    :Returns:
        Return the arguments of :func:`gen_adel_input_data`. 
    
    :Returns Type:
        tuple
    
    '''
    import imp

    inputs = imp.load_source('inputs', inputs_filepath)

    dynT_user = pandas.read_csv(inputs.dynT_user)
    dimT_user = pandas.read_csv(inputs.dimT_user)
    plants_number = inputs.plants_number
    plants_density = inputs.plants_density
    decide_child_axis_probabilities = inputs.decide_child_axis_probabilities 
    MS_leaves_number_probabilities = inputs.MS_leaves_number_probabilities
    ears_density = inputs.ears_density
    GL_number = inputs.GL_number
    delais_TT_stop_del_axis = inputs.delais_TT_stop_del_axis
    TT_col_break = inputs.TT_col_break
    try:
        inner_params = inputs.inner_params
    except:
        inner_params = {}
    
    return (dynT_user, dimT_user, plants_number, plants_density, decide_child_axis_probabilities, 
            MS_leaves_number_probabilities, ears_density, GL_number, 
            delais_TT_stop_del_axis, TT_col_break, inner_params)
