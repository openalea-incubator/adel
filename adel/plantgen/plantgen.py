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

from adel.plantgen import axeT, dimT, dynT, phenT, tools, params

class DataCompleteness:
    '''
    This enumerate defines the different degrees of completeness that the data 
    documented by the user can have.
    
    .. seealso:: :func:`gen_adel_input_data`  
    '''
    MIN=1
    SHORT=2
    FULL=3


def gen_adel_input_data_from_min(dynT_user={'a_cohort': 0.0102, 'TT_col_0': -0.771289027, 'n0': 4.871559739, 'n1': 3.24283148, 'n2': 5.8},
                                 TT_col_nff={'1': 1078.0, '4': 1148.0, '5': 1158.0, '6': 1168.0, '7': 1178.0},
                                 dimT_user=None,
                                 plant_number=100, 
                                 decide_child_cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.983, '6': 0.817, '7': 0.117}, 
                                 MS_leaves_number_probabilities={'10': 0.145, '11': 0.818, '12': 0.037, '13': 0.0, '14': 0.0},
                                 TT_bolting=500.0,
                                 final_axes_density=250,
                                 GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                                 delais_TT_stop_del_axis=600,
                                 TT_col_break=0.0):
    '''
    Generate ADEL input data from a *MIN* set of data. This a convenience 
    function to be used from VisuAlea. 
    
    The *MIN* set of data is represented by *dynT_user*, *TT_col_nff* and *dimT_user*. 
    See :ref:`plantgen` for an example of how to set these parameters properly.
    
    :Parameters:
    
        - `dynT_user` (:class:`dict`) - the leaf dynamic parameters set by the user. See 
          :ref:`dynT_user_MIN <dynT_user_MIN>`
                      
          *dynT_user* must be a dict with the following shape:: 
        
              {'a_cohort': a_cohort, 'TT_col_0': TT_col_0, 
               'n0': n0, 'n1': n1, 'n2': n2}
        
          where ``a_cohort``, ``TT_col_0``, ``n0``, ``n1`` and ``n2`` are floats.
                
        - `TT_col_nff` (:class:`dict`) - the thermal time when Haun Stage is equal to 
          *Nff*, for the main stem.
          
          *TT_col_nff* must be a dict with the following shape:: 
        
              {'1': value_1, '4': value_4, '5': value_5, 
               '6': value_6, '7': value_7}
        
          where ``value_\*`` are floats. 
        
        - `dimT_user` (:class:`pandas.DataFrame`) - the dimensions of the organs set by 
          the user. See :ref:`dimT_user_MIN <dimT_user_MIN>`.
          
          *dimT_user* must be a pandas.Dataframe with the 
          following columns: *index_phytomer*, *L_blade*, *W_blade*, *L_sheath*, 
          *W_sheath*, *L_internode*, *W_internode*.
          The values can be either integers or floats.
              
        - `plant_number` (:class:`int`) - the number of plants to be generated.
        
        - `decide_child_cohort_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          for each child cohort the probability of emergence of an axis when the parent 
          axis is present. The keys are the numbers of the child cohorts and the 
          values are the probabilities.
        
        - `MS_leaves_number_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          the probability distribution of the final number of main stem leaves. 
          The keys are the final numbers of main stem leaves, and the values are 
          the probability distribution.
          
        - `TT_bolting` (:class:`int`) - date in thermal time at which the bolting starts.
        
        - `final_axes_density` (:class:`int`) - the final number of axes which have an 
          ear, per square meter.
        
        - `GL_number` (:class:`dict` of :class:`float`::class:`float`) - the GL decimal numbers measured at 
          several thermal times (including the senescence end). The keys are the 
          thermal times, and the values are the GL decimal numbers.
        
        - `delais_TT_stop_del_axis` (:class:`int`) - This variable represents the time in 
          thermal time between an axis stop growing and its disappearance (it 
          concerns only the axes that do not regress and which do not produce any 
          cob).
        
        - `TT_col_break` (:class:`float`) - the thermal time when the rate of Haun Stage 
          is changing. If phyllochron is constant, then *TT_col_break* is null.
        
    :Returns:
        Return the following dataframes: :ref:`axeT <axeT>`, :ref:`dimT <dimT>`, 
        :ref:`phenT <phenT>`, :ref:`phenT_abs <phenT_abs>`, :ref:`dimT_abs <dimT_abs>`, 
        :ref:`dynT <dynT>`, :ref:`phenT_first <phenT_first>`, 
        :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`, :ref:`tilleringT <tilleringT>`, 
        :ref:`cohortT <cohortT>`
    
    :Returns Type:
        tuple of :class:`pandas.DataFrame`
        
    .. seealso:: :ref:`plantgen`
                 :func:`gen_adel_input_data`
                 :mod:`alinea.adel.plantgen.axeT`
                 :mod:`alinea.adel.plantgen.dimT`
                 :mod:`alinea.adel.plantgen.dynT`
                 :mod:`alinea.adel.plantgen.params`
                 :mod:`alinea.adel.plantgen.phenT`
                 :mod:`alinea.adel.plantgen.tools`
                 
    .. warning:: the type of the arguments is checked as follows:

             .. list-table::
                 :widths: 10 50
                 :header-rows: 1
            
                 * - Argument
                   - Type
                 * - *dynT_user* 
                   - :class:`dict`
                 * - *TT_col_nff* 
                   - :class:`dict`  
                 * - *dimT_user* 
                   - :class:`pandas.DataFrame`
                 * - *plant_number* 
                   - :class:`int`
                 * - *decide_child_cohort_probabilities* 
                   - :class:`dict`
                 * - *MS_leaves_number_probabilities* 
                   - :class:`dict`
                 * - *TT_bolting* 
                   - :class:`float`
                 * - *final_axes_density* 
                   - :class:`int`
                 * - *GL_number* 
                   - :class:`dict`
                 * - *delais_TT_stop_del_axis* 
                   - :class:`int`
                 * - *TT_col_break* 
                   - :class:`float`

    '''    
    tools.checkValidity(isinstance(dynT_user, dict))
    tools.checkValidity(isinstance(TT_col_nff, dict))
    tools.checkValidity(isinstance(dimT_user, pandas.DataFrame))
    tools.checkValidity(isinstance(plant_number, int))
    tools.checkValidity(isinstance(decide_child_cohort_probabilities, dict))
    tools.checkValidity(isinstance(MS_leaves_number_probabilities, dict))
    tools.checkValidity(isinstance(TT_bolting, float))
    tools.checkValidity(isinstance(final_axes_density, int))
    tools.checkValidity(isinstance(GL_number, dict))
    tools.checkValidity(isinstance(delais_TT_stop_del_axis, int))
    tools.checkValidity(isinstance(TT_col_break, float))
            
    # check dynT_user validity
    expected_dynT_user_keys_value_types = {'a_cohort': float, 
                                           'TT_col_0': float, 
                                           'n0': float,
                                           'n1': float,
                                           'n2': float}
    dynT_user_keys_value_types = dict(zip(dynT_user.keys(), 
                                          [type(value) for value in dynT_user.values()]))
    tools.checkValidity(expected_dynT_user_keys_value_types == dynT_user_keys_value_types)
    # check TT_col_nff validity
    expected_TT_col_nff_keys_value_types = {'1': float, '4': float, '5': float, '6': float, '7': float}
    TT_col_nff_keys_value_types = dict(zip(TT_col_nff.keys(), 
                                          [type(value) for value in TT_col_nff.values()]))
    tools.checkValidity(expected_TT_col_nff_keys_value_types == TT_col_nff_keys_value_types)
    # check dimT_user validity
    tools.checkValidity(isinstance(dimT_user, pandas.DataFrame))
    expected_dimT_user_columns = ['index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
    tools.checkValidity(dimT_user.columns.tolist() == expected_dimT_user_columns)
    tools.checkValidity(dimT_user.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all())
    tools.checkValidity(dimT_user['index_phytomer'].unique().size == dimT_user['index_phytomer'].size)
    dynT_user = dynT_user.copy()
    dynT_user['TT_col_nff'] = TT_col_nff
    return gen_adel_input_data(dynT_user, dimT_user, plant_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities, TT_bolting, final_axes_density, GL_number, delais_TT_stop_del_axis, TT_col_break, DataCompleteness.MIN, DataCompleteness.MIN)


def gen_adel_input_data_from_short(dynT_user,
                                    dimT_user,
                                    plant_number=100, 
                                    decide_child_cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.983, '6': 0.817, '7': 0.117}, 
                                    MS_leaves_number_probabilities={'10': 0.145, '11': 0.818, '12': 0.037, '13': 0.0, '14': 0.0},
                                    TT_bolting=500.0,
                                    final_axes_density=250,
                                    GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                                    delais_TT_stop_del_axis=600,
                                    TT_col_break=0.0):
    '''
    Generate ADEL input data from a *SHORT* set of data. This a convenience 
    function to be used from VisuAlea. 
    
    The *SHORT* set of data is represented by *dynT_user* and *dimT_user*.
    See :ref:`plantgen` for an example of how to set these parameters properly.
    
    :Parameters:
    
        - `dynT_user` (:class:`pandas.DataFrame`) - the leaf dynamic parameters 
          set by the user. See :ref:`dynT_user_SHORT <dynT_user_SHORT>`.
               
          *dynT_user* must be a pandas.Dataframe with the 
          following columns: *N_cohort*, *a_cohort*, *TT_col_0*, *TT_col_nff*, *n0*, *n1*, *n2*.
          The values can be either integers or floats.
        
        - `dimT_user` (:class:`pandas.DataFrame`) - the dimensions of the organs set by 
          the user. See :ref:`dimT_user_SHORT <dimT_user_SHORT>`
          
          *dimT_user* must be a pandas.Dataframe with the 
          following columns: *id_axis*, *index_phytomer*, *L_blade*, *W_blade*, 
          *L_sheath*, *W_sheath*, *L_internode*, *W_internode*.
          The values can be either integers or floats.
              
        - `plant_number` (:class:`int`) - the number of plants to be generated.
        
        - `decide_child_cohort_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          for each child cohort the probability of emergence of an axis when the parent 
          axis is present. The keys are the numbers of the child cohorts and the 
          values are the probabilities.
        
        - `MS_leaves_number_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          the probability distribution of the final number of main stem leaves. 
          The keys are the final numbers of main stem leaves, and the values are 
          the probability distribution.
          
        - `TT_bolting` (:class:`int`) - date in thermal time at which the bolting starts.
        
        - `final_axes_density` (:class:`int`) - the final number of axes which have an 
          ear, per square meter.
        
        - `GL_number` (:class:`dict` of :class:`float`::class:`float`) - the GL decimal numbers measured at 
          several thermal times (including the senescence end). The keys are the 
          thermal times, and the values are the GL decimal numbers.
        
        - `delais_TT_stop_del_axis` (:class:`int`) - This variable represents the time in 
          thermal time between an axis stop growing and its disappearance (it 
          concerns only the axes that do not regress and which do not produce any 
          cob).
        
        - `TT_col_break` (:class:`float`) - the thermal time when the rate of Haun Stage 
          is changing. If phyllochron is constant, then *TT_col_break* is null.
        
    :Returns:
        Return the following dataframes: :ref:`axeT <axeT>`, :ref:`dimT <dimT>`, 
        :ref:`phenT <phenT>`, :ref:`phenT_abs <phenT_abs>`, :ref:`dimT_abs <dimT_abs>`, 
        :ref:`dynT <dynT>`, :ref:`phenT_first <phenT_first>`, 
        :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`, :ref:`tilleringT <tilleringT>`,
        :ref:`cohortT <cohortT>`
    
    :Returns Type:
        tuple of :class:`pandas.DataFrame`
        
    .. seealso:: :ref:`plantgen`
                 :func:`gen_adel_input_data`
                 :mod:`alinea.adel.plantgen.axeT`
                 :mod:`alinea.adel.plantgen.dimT`
                 :mod:`alinea.adel.plantgen.dynT`
                 :mod:`alinea.adel.plantgen.params`
                 :mod:`alinea.adel.plantgen.phenT`
                 :mod:`alinea.adel.plantgen.tools`
                 
    .. warning:: the arguments are checked as follows:

         .. list-table::
             :widths: 10 10 30
             :header-rows: 1
        
             * - Argument
               - Type
               - Constraint
             * - *dynT_user* 
               - :class:`pandas.DataFrame`
               - **MUST** contain a row of data for each possible cohort. 
                 See *decide_child_cohort_probabilities*.
             * - *dimT_user* 
               - :class:`pandas.DataFrame`
               - **MUST** contain a row of data for each possible cohort. 
                 See *decide_child_cohort_probabilities*.
             * - *plant_number* 
               - :class:`int`
               - None
             * - *decide_child_cohort_probabilities* 
               - :class:`dict`
               - None
             * - *MS_leaves_number_probabilities* 
               - :class:`dict`
               - None
             * - *TT_bolting* 
               - :class:`float`
               - None
             * - *final_axes_density* 
               - :class:`int`
               - None
             * - *GL_number* 
               - :class:`dict`
               - None
             * - *delais_TT_stop_del_axis* 
               - :class:`int`
               - None
             * - *TT_col_break* 
               - :class:`float`
               - None

    '''    
    tools.checkValidity(isinstance(dynT_user, pandas.DataFrame))
    tools.checkValidity(isinstance(dimT_user, pandas.DataFrame))
    tools.checkValidity(isinstance(plant_number, int))
    tools.checkValidity(isinstance(decide_child_cohort_probabilities, dict))
    tools.checkValidity(isinstance(MS_leaves_number_probabilities, dict))
    tools.checkValidity(isinstance(TT_bolting, float))
    tools.checkValidity(isinstance(final_axes_density, int))
    tools.checkValidity(isinstance(GL_number, dict))
    tools.checkValidity(isinstance(delais_TT_stop_del_axis, int))
    tools.checkValidity(isinstance(TT_col_break, float))

    return gen_adel_input_data(dynT_user, dimT_user, plant_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities, TT_bolting, final_axes_density, GL_number, delais_TT_stop_del_axis, TT_col_break, DataCompleteness.SHORT, DataCompleteness.SHORT)
    

def gen_adel_input_data_from_full(dynT_user,
                                    dimT_user,
                                    plant_number=100, 
                                    decide_child_cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.983, '6': 0.817, '7': 0.117}, 
                                    MS_leaves_number_probabilities={'10': 0.145, '11': 0.818, '12': 0.037, '13': 0.0, '14': 0.0},
                                    TT_bolting=500.0,
                                    final_axes_density=250,
                                    GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                                    delais_TT_stop_del_axis=600,
                                    TT_col_break=0.0):
    '''
    Generate ADEL input data from a *FULL* set of data. This a convenience 
    function to be used from VisuAlea. 
    
    *FULL* set of data is represented by *dynT_user* and *dimT_user*.
    See :ref:`plantgen` for an example of how to set these parameters properly.
    
    :Parameters:
    
        - `dynT_user` (:class:`pandas.DataFrame`) - the leaf dynamic parameters 
          set by the user. See :ref:`dynT_user_FULL <dynT_user_FULL>`.
               
          *dynT_user* must be a pandas.Dataframe with the 
          following columns: *N_cohort*, *Nff*, *a_cohort*, *TT_col_0*, *TT_col_nff*, *n0*, *n1*, *n2*.
          The values can be either integers or floats.
        
        - `dimT_user` (:class:`pandas.DataFrame`) - the dimensions of the organs set by 
          the user. See :ref:`dimT_user_FULL <dimT_user_FULL>`.
          
          *dimT_user* must be a pandas.Dataframe with the 
          following columns: *id_dim*, *index_phytomer*, *L_blade*, *W_blade*, 
          *L_sheath*, *W_sheath*, *L_internode*, *W_internode*.
          The values can be either integers or floats.
              
        - `plant_number` (:class:`int`) - the number of plants to be generated.
        
        - `decide_child_cohort_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          for each child cohort the probability of emergence of an axis when the parent 
          axis is present. The keys are the numbers of the child cohorts and the 
          values are the probabilities.
        
        - `MS_leaves_number_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          the probability distribution of the final number of main stem leaves. 
          The keys are the final numbers of main stem leaves, and the values are 
          the probability distribution.
          
        - `TT_bolting` (:class:`int`) - date in thermal time at which the bolting starts.
        
        - `final_axes_density` (:class:`int`) - the final number of axes which have an 
          ear, per square meter.
        
        - `GL_number` (:class:`dict` of :class:`float`::class:`float`) - the GL decimal numbers measured at 
          several thermal times (including the senescence end). The keys are the 
          thermal times, and the values are the GL decimal numbers.
        
        - `delais_TT_stop_del_axis` (:class:`int`) - This variable represents the time in 
          thermal time between an axis stop growing and its disappearance (it 
          concerns only the axes that do not regress and which do not produce any 
          cob).
        
        - `TT_col_break` (:class:`float`) - the thermal time when the rate of Haun Stage 
          is changing. If phyllochron is constant, then *TT_col_break* is null.
        
    :Returns:
        Return the following dataframes: :ref:`axeT <axeT>`, :ref:`dimT <dimT>`, 
        :ref:`phenT <phenT>`, :ref:`phenT_abs <phenT_abs>`, :ref:`dimT_abs <dimT_abs>`, 
        :ref:`dynT <dynT>`, :ref:`phenT_first <phenT_first>`, 
        :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`, :ref:`tilleringT <tilleringT>`,
        :ref:`cohortT <cohortT>`.
    
    :Returns Type:
        tuple of :class:`pandas.DataFrame`
        
    .. seealso:: :ref:`plantgen`
                 :func:`gen_adel_input_data`
                 :mod:`alinea.adel.plantgen.axeT`
                 :mod:`alinea.adel.plantgen.dimT`
                 :mod:`alinea.adel.plantgen.dynT`
                 :mod:`alinea.adel.plantgen.params`
                 :mod:`alinea.adel.plantgen.phenT`
                 :mod:`alinea.adel.plantgen.tools`
                 
    .. warning:: the arguments are checked as follows:

         .. list-table::
             :widths: 10 10 30
             :header-rows: 1
        
             * - Argument
               - Type
               - Constraint
             * - *dynT_user*
               - :class:`pandas.DataFrame`
               - **MUST** contain a row of data for each possible phytomer of the 
                 most frequent axis of the main stem, and for each possible cohort. 
                 See *decide_child_cohort_probabilities* and *MS_leaves_number_probabilities*.
             * - *dimT_user*
               - :class:`pandas.DataFrame`
               - **MUST** contain a row of data for each possible phytomer of the 
                 most frequent axis of the main stem, and for each possible cohort. 
                 See *decide_child_cohort_probabilities* and *MS_leaves_number_probabilities*.
             * - *plant_number* 
               - :class:`int`
               - None
             * - *decide_child_cohort_probabilities* 
               - :class:`dict`
               - None
             * - *MS_leaves_number_probabilities* 
               - :class:`dict`
               - None
             * - *TT_bolting* 
               - :class:`float`
               - None
             * - *final_axes_density* 
               - :class:`int`
               - None
             * - *GL_number* 
               - :class:`dict`
               - None
             * - *delais_TT_stop_del_axis* 
               - :class:`int`
               - None
             * - *TT_col_break* 
               - :class:`float`
               - None    
                      
    '''    
    tools.checkValidity(isinstance(dynT_user, pandas.DataFrame))
    tools.checkValidity(isinstance(dimT_user, pandas.DataFrame))
    tools.checkValidity(isinstance(plant_number, int))
    tools.checkValidity(isinstance(decide_child_cohort_probabilities, dict))
    tools.checkValidity(isinstance(MS_leaves_number_probabilities, dict))
    tools.checkValidity(isinstance(TT_bolting, float))
    tools.checkValidity(isinstance(final_axes_density, int))
    tools.checkValidity(isinstance(GL_number, dict))
    tools.checkValidity(isinstance(delais_TT_stop_del_axis, int))
    tools.checkValidity(isinstance(TT_col_break, float))

    return gen_adel_input_data(dynT_user, dimT_user, plant_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities, TT_bolting, final_axes_density, GL_number, delais_TT_stop_del_axis, TT_col_break, DataCompleteness.FULL, DataCompleteness.FULL)


def gen_adel_input_data(dynT_user,
                        dimT_user,
                        plant_number=100, 
                        decide_child_cohort_probabilities={'3': 0.0, '4': 0.900, '5': 0.983, '6': 0.817, '7': 0.117}, 
                        MS_leaves_number_probabilities={'10': 0.145, '11': 0.818, '12': 0.037, '13': 0.0, '14': 0.0},
                        TT_bolting=500.0,
                        final_axes_density=250,
                        GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                        delais_TT_stop_del_axis=600,
                        TT_col_break=0.0,
                        dynT_user_completeness=DataCompleteness.MIN,
                        dimT_user_completeness=DataCompleteness.MIN                        
                        ):
    '''
    Create the dataframes which contain the plant data to be used as input for 
    generating plot with ADEL, and some other dataframes for debugging purpose.
    See :ref:`adel_input` for a description of the input tables expected by ADEL, 
    and :ref:`plantgen` for a description of the dataframes created for debug. 
    
    Different degrees of completeness of data provided by the user are acceptable. 
    The user must specify the degree of completeness selecting a value within the 
    enumerate :class:`DataCompleteness`.
    
    The dataframes are created as follows:
        * initialization of the following dataframes:
            * *axeT_tmp*, calling :func:`alinea.adel.plantgen.axeT.create_axeT_tmp`,
            * *dynT_tmp*, calling :func:`alinea.adel.plantgen.dynT.create_dynT_tmp`.
            * *dimT_tmp*, calling :func:`alinea.adel.plantgen.dimT.create_dimT_tmp`
            * :ref:`tilleringT <tilleringT>`, calling :func:`alinea.adel.plantgen.axeT.create_tilleringT`
            * :ref:`cohortT <cohortT>`, calling :func:`alinea.adel.plantgen.axeT.create_cohortT`
        * filling of the dataframes set by the user:
            * *dynT_user*, according to *dynT_user_completeness*
            * *dimT_user*, according to *dimT_user_completeness*
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
            * the :ref:`cohortT <cohortT>`,
            * the :ref:`phenT_abs <phenT_abs>`,
            * the :ref:`dimT_abs <dimT_abs>`,
            * the :ref:`dynT <dynT>`, 
            * the :ref:`phenT_first <phenT_first>`,
            * the :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`,
        
    :Parameters:
    
        - `dynT_user` (:class:`dict` | :class:`pandas.DataFrame`) - the leaf dynamic 
          parameters set by the user.
          The type and the content depend on the *dynT_user_completeness*.
          See :ref:`levels_of_completeness`.       
                
        - `dimT_user` (:class:`pandas.DataFrame`) - the dimensions of the organs 
          set by the user. 
          The content depends on the *dimT_user_completeness* argument. 
          See :ref:`levels_of_completeness`.
              
        - `plant_number` (:class:`int`) - the number of plants to be generated.
        
        - `decide_child_cohort_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          for each child cohort the probability of emergence of an axis when the parent 
          axis is present. The keys are the numbers of the child cohorts and the 
          values are the probabilities.
        
        - `MS_leaves_number_probabilities` (:class:`dict` of :class:`str`::class:`float`) - 
          the probability distribution of the final number of main stem leaves. 
          The keys are the final numbers of main stem leaves, and the values are 
          the probabilities distribution.
          
        - `TT_bolting` (:class:`int`) - date in thermal time at which the bolting starts.
        
        - `final_axes_density` (:class:`int`) - the final number of axes which have an 
          ear, per square meter.
        
        - `GL_number` (:class:`dict` of :class:`float`::class:`float`) - the GL decimal numbers measured at 
          several thermal times (including the senescence end). The keys are the 
          thermal times, and the values are the GL decimal numbers.
        
        - `delais_TT_stop_del_axis` (:class:`int`) - This variable represents the time in 
          thermal time between an axis stop growing and its disappearance (it 
          concerns only the axes that do not regress and which do not produce any 
          cob).
        
        - `TT_col_break` (:class:`float`) - the thermal time when the rate of Haun Stage 
          is changing. If phyllochron is constant, then *TT_col_break* is null.
        
        - `dynT_user_completeness` (:class:`DataCompleteness`) - the level of 
          completeness of the *dynT_user* set by the user. 
        
        - `dimT_user_completeness` (:class:`DataCompleteness`) - the level of completeness of the 
          *dimT_user* set by the user. 
        
    :Returns:
        Return the following dataframes: :ref:`axeT <axeT>`, :ref:`dimT <dimT>`, 
        :ref:`phenT <phenT>`, :ref:`phenT_abs <phenT_abs>`, :ref:`dimT_abs <dimT_abs>`, 
        :ref:`dynT <dynT>`, :ref:`phenT_first <phenT_first>`, :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`, 
        :ref:`tilleringT <tilleringT>`, :ref:`cohortT <cohortT>`.
    
    :Returns Type:
        tuple of :class:`pandas.DataFrame`
        
    .. seealso:: :class:`DataCompleteness`
                 :mod:`alinea.adel.plantgen.axeT`
                 :mod:`alinea.adel.plantgen.dimT`
                 :mod:`alinea.adel.plantgen.dynT`
                 :mod:`alinea.adel.plantgen.params`
                 :mod:`alinea.adel.plantgen.phenT`
                 :mod:`alinea.adel.plantgen.tools`
                 
    .. warning:: the type of the arguments is checked as follows:
    
                 .. list-table::
                     :widths: 10 20 20 20
                     :header-rows: 1
                
                     * - Argument
                       - MIN
                       - SHORT
                       - FULL
                     * - *dynT_user*
                       - :class:`dict`
                       - :class:`pandas.DataFrame` which contains a row of data 
                         for each possible cohort. See *decide_child_cohort_probabilities*.
                       - :class:`pandas.DataFrame` which contains a row of data 
                         for each possible leaves number of the most frequent axis 
                         of the main stem, and for each possible cohort. 
                         See *decide_child_cohort_probabilities* and *MS_leaves_number_probabilities*.
                     * - *dimT_user*
                       - :class:`pandas.DataFrame`
                       - :class:`pandas.DataFrame` which contains a row of data 
                         for each possible cohort. See *decide_child_cohort_probabilities*.
                       - :class:`pandas.DataFrame` which contains a row of data 
                         for each possible leaves number of the most frequent axis 
                         of the main stem, and for each possible cohort.
                         See *decide_child_cohort_probabilities* and *MS_leaves_number_probabilities*. 
                     * - *plant_number* 
                       - :class:`int`
                       - :class:`int`
                       - :class:`int`
                     * - *decide_child_cohort_probabilities* 
                       - :class:`dict`
                       - :class:`dict`
                       - :class:`dict`
                     * - *MS_leaves_number_probabilities* 
                       - :class:`dict`
                       - :class:`dict`
                       - :class:`dict`
                     * - *TT_bolting* 
                       - :class:`float`
                       - :class:`float`
                       - :class:`float`
                     * - *final_axes_density* 
                       - :class:`int`
                       - :class:`int`
                       - :class:`int`
                     * - *GL_number* 
                       - :class:`dict`
                       - :class:`dict`
                       - :class:`dict`
                     * - *delais_TT_stop_del_axis* 
                       - :class:`int`
                       - :class:`int`
                       - :class:`int`
                     * - *TT_col_break* 
                       - :class:`float`
                       - :class:`float`
                       - :class:`float`
                     * - *dynT_user_completeness* 
                       - :attr:`DataCompleteness.MIN`
                       - :attr:`DataCompleteness.SHORT`
                       - :attr:`DataCompleteness.FULL`
                     * - *dimT_user_completeness* 
                       - :attr:`DataCompleteness.MIN`
                       - :attr:`DataCompleteness.SHORT`
                       - :attr:`DataCompleteness.FULL`
    
    '''
    tools.checkValidity(isinstance(dynT_user, (dict, pandas.DataFrame)))
    tools.checkValidity(isinstance(dimT_user, pandas.DataFrame))
    tools.checkValidity(isinstance(plant_number, int))
    tools.checkValidity(isinstance(decide_child_cohort_probabilities, dict))
    tools.checkValidity(isinstance(MS_leaves_number_probabilities, dict))
    tools.checkValidity(isinstance(TT_bolting, float))
    tools.checkValidity(isinstance(final_axes_density, int))
    tools.checkValidity(isinstance(GL_number, dict))
    tools.checkValidity(isinstance(delais_TT_stop_del_axis, int))
    tools.checkValidity(isinstance(TT_col_break, float))
    tools.checkValidity(dynT_user_completeness in DataCompleteness.__dict__.values())
    tools.checkValidity(dimT_user_completeness in DataCompleteness.__dict__.values())
    
    possible_cohorts = \
        set([idx_of_cohort for (idx_of_cohort, probability) in
             decide_child_cohort_probabilities.iteritems() if probability != 0.0])
        
    # check plant_number, decide_child_cohort_probabilities and final_axes_density validity
    theoretical_cohorts_cardinalities = tools.calculate_theoretical_cohorts_cardinalities(plant_number, 
                                                                                          decide_child_cohort_probabilities,
                                                                                          params.FIRST_CHILD_DELAY)
    theoretical_cardinalities_sum = sum(theoretical_cohorts_cardinalities.values())
    tools.checkValidity(final_axes_density <= int(theoretical_cardinalities_sum))
    
    # check dynT_user validity
    if dynT_user_completeness == DataCompleteness.MIN:
        expected_dynT_user_keys_value_types = {'a_cohort': float, 
                                                 'TT_col_0': float, 
                                                 'TT_col_nff': dict, 
                                                 'n0': float,
                                                 'n1': float,
                                                 'n2': float}
        dynT_user_keys_value_types = dict(zip(dynT_user.keys(), 
                                              [type(value) for value in dynT_user.values()]))
        tools.checkValidity(expected_dynT_user_keys_value_types == dynT_user_keys_value_types)
        # check the validity of the dict associated to 'TT_col_nff'
        expected_TT_col_nff_keys_value_types = {'1': float, '4': float, '5': float, '6': float, '7': float}
        TT_col_nff_keys_value_types = dict(zip(dynT_user['TT_col_nff'].keys(), 
                                               [type(value) for value in dynT_user['TT_col_nff'].values()]))
        tools.checkValidity(expected_TT_col_nff_keys_value_types == TT_col_nff_keys_value_types)
    elif dynT_user_completeness == DataCompleteness.SHORT:
        tools.checkValidity(isinstance(dynT_user, pandas.DataFrame))
        expected_dynT_user_columns = ['N_cohort', 'a_cohort', 'TT_col_0', 'TT_col_nff', 'n0', 'n1', 'n2']
        tools.checkValidity(dynT_user.columns.tolist() == expected_dynT_user_columns)
        tools.checkValidity(dynT_user.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all())
        tools.checkValidity(dynT_user['N_cohort'].unique().size == dynT_user['N_cohort'].size)
        available_cohorts = set(dynT_user['N_cohort'].astype(int).astype(str).tolist())
        tools.checkValidity(possible_cohorts.issubset(available_cohorts))
    elif dynT_user_completeness == DataCompleteness.FULL:
        tools.checkValidity(isinstance(dynT_user, pandas.DataFrame))
        expected_dynT_user_columns = ['N_cohort', 'Nff', 'a_cohort', 'TT_col_0', 'TT_col_nff', 'n0', 'n1', 'n2']
        tools.checkValidity(dynT_user.columns.tolist() == expected_dynT_user_columns)
        tools.checkValidity(dynT_user.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all())
        grouped = dynT_user.groupby(['N_cohort', 'Nff'])
        tools.checkValidity(len(grouped.groups) == dynT_user.index.size   )
        available_cohorts = set(dynT_user['N_cohort'].astype(int).astype(str).tolist())
        tools.checkValidity(possible_cohorts.issubset(available_cohorts) )
        
    # check dimT_user validity
    if dimT_user_completeness == DataCompleteness.MIN:
        tools.checkValidity(isinstance(dimT_user, pandas.DataFrame))
        expected_dimT_user_columns = ['index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        tools.checkValidity(dimT_user.columns.tolist() == expected_dimT_user_columns)
        tools.checkValidity(dimT_user.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all())
        tools.checkValidity(dimT_user['index_phytomer'].unique().size == dimT_user['index_phytomer'].size)
    elif dimT_user_completeness == DataCompleteness.SHORT:
        tools.checkValidity(isinstance(dimT_user, pandas.DataFrame))
        expected_dimT_user_columns = ['id_axis', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        tools.checkValidity(dimT_user.columns.tolist() == expected_dimT_user_columns)
        tools.checkValidity(dimT_user.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all())
        grouped = dimT_user.groupby(['id_axis', 'index_phytomer'])
        tools.checkValidity(len(grouped.groups) == dimT_user.index.size    )
        available_cohorts = set(dimT_user['id_axis'].astype(int).astype(str).tolist())
        tools.checkValidity(possible_cohorts.issubset(available_cohorts))
    elif dimT_user_completeness == DataCompleteness.FULL:
        tools.checkValidity(isinstance(dimT_user, pandas.DataFrame))
        expected_dimT_user_columns = ['id_dim', 'index_phytomer', 'L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
        tools.checkValidity(dimT_user.columns.tolist() == expected_dimT_user_columns)
        tools.checkValidity(dimT_user.dtypes.isin([np.dtype(np.int64), np.dtype(np.float64)]).all())
        grouped = dimT_user.groupby(['id_dim', 'index_phytomer'])
        tools.checkValidity(len(grouped.groups) == dimT_user.index.size)
        available_cohorts = set([str(id_dim)[:-2] for id_dim in dimT_user['id_dim'].astype(int).tolist()])
        tools.checkValidity(possible_cohorts.issubset(available_cohorts))
    
    # 2. first step of the fit process
    if dynT_user_completeness == DataCompleteness.MIN:
        TT_flag_leaf_ligulation = dynT_user['TT_col_nff']['1']
    else:
        TT_flag_leaf_ligulation = dynT_user['TT_col_nff'][dynT_user.first_valid_index()]
    
    (axeT_tmp_dataframe, 
     dimT_tmp_dataframe, 
     dynT_tmp_dataframe, 
     tilleringT_dataframe,
     cohortT_dataframe) = _gen_adel_input_data_first(plant_number, 
                                                     decide_child_cohort_probabilities, 
                                                     MS_leaves_number_probabilities, 
                                                     TT_bolting, 
                                                     TT_flag_leaf_ligulation, 
                                                     final_axes_density,
                                                     theoretical_cohorts_cardinalities) 
    
    # 3. complete dynT_user
    if dynT_user_completeness == DataCompleteness.MIN:
        dynT_tmp_dataframe.ix[0]['TT_col_break'] = TT_col_break
        dynT_tmp_dataframe.ix[0]['a_cohort'] = dynT_user['a_cohort']
        dynT_tmp_dataframe.ix[0]['TT_col_0'] = dynT_user['TT_col_0']
        TT_col_nff_keys = dynT_user['TT_col_nff'].keys()
        TT_col_nff_keys.sort()
        first_TT_col_nff = dynT_user['TT_col_nff'][TT_col_nff_keys[0]]
        for N_cohort, dynT_tmp_dataframe_grouped_by_N_cohort in dynT_tmp_dataframe.groupby('N_cohort'):
            N_cohort_int = int(N_cohort)
            current_TT_col_nff = dynT_user['TT_col_nff'][str(N_cohort_int)]
            current_dTT_MS_cohort = current_TT_col_nff - first_TT_col_nff
            index_to_set = dynT_tmp_dataframe_grouped_by_N_cohort.index[0]
            dynT_tmp_dataframe['TT_col_nff'][index_to_set] = current_TT_col_nff
            dynT_tmp_dataframe['dTT_MS_cohort'][index_to_set] = current_dTT_MS_cohort
        dynT_tmp_dataframe.ix[0]['n0'] = dynT_user['n0']
        dynT_tmp_dataframe.ix[0]['n1'] = dynT_user['n1']
        dynT_tmp_dataframe.ix[0]['n2'] = dynT_user['n2']
    elif dynT_user_completeness == DataCompleteness.SHORT:
        user_grouped = dynT_user.groupby('N_cohort')
        for N_cohort, generated_group in dynT_tmp_dataframe.groupby('N_cohort'):
            if N_cohort not in user_grouped.groups:
                continue
            user_group = user_grouped.get_group(N_cohort)
            index_to_get = user_group.index[0]
            index_to_set = generated_group.index[0]
            if N_cohort == 1.0:
                first_TT_col_nff = user_group['TT_col_nff'][index_to_get]
            columns_to_set = dynT_user.columns
            dynT_tmp_dataframe.ix[index_to_set][columns_to_set] = dynT_user.ix[index_to_get]
            dynT_tmp_dataframe.ix[index_to_set]['TT_col_break'] = TT_col_break
            current_TT_col_nff = user_group['TT_col_nff'][index_to_get]
            current_dTT_MS_cohort = current_TT_col_nff - first_TT_col_nff
            dynT_tmp_dataframe.ix[index_to_set]['dTT_MS_cohort'] = current_dTT_MS_cohort
    elif dynT_user_completeness == DataCompleteness.FULL:
        user_grouped_N_cohort = dynT_user.groupby('N_cohort')
        for N_cohort, N_cohort_generated_group in dynT_tmp_dataframe.groupby('N_cohort'):
            if N_cohort not in user_grouped_N_cohort.groups:
                continue
            user_N_cohort_group = user_grouped_N_cohort.get_group(N_cohort)
            i = 0
            most_frequent_Nff = N_cohort_generated_group['Nff'][N_cohort_generated_group.index[i]]
            while user_N_cohort_group[user_N_cohort_group['Nff'] == most_frequent_Nff].index.size == 0:
                i += 1
                most_frequent_Nff = N_cohort_generated_group['Nff'][N_cohort_generated_group.index[i]]
            N_cohort_index_to_get = user_N_cohort_group[user_N_cohort_group['Nff'] == most_frequent_Nff].index[0]
            N_cohort_index_to_set = N_cohort_generated_group[N_cohort_generated_group['Nff'] == most_frequent_Nff].index[0]
            if N_cohort == 1.0:
                first_TT_col_nff = user_N_cohort_group['TT_col_nff'][N_cohort_index_to_get]
                most_frequent_MS_a_cohort = dynT_user['a_cohort'][N_cohort_index_to_get]
            current_TT_col_nff = user_N_cohort_group['TT_col_nff'][N_cohort_index_to_get]
            most_frequent_axis_dTT_MS_cohort = current_TT_col_nff - first_TT_col_nff
            dynT_tmp_dataframe.ix[N_cohort_index_to_set]['dTT_MS_cohort'] = most_frequent_axis_dTT_MS_cohort
            highest_cardinality = N_cohort_generated_group['cardinality'][N_cohort_index_to_set]
            user_grouped_Nff = user_N_cohort_group.groupby('Nff')
            for Nff, Nff_generated_group in N_cohort_generated_group.groupby('Nff'):
                if Nff not in user_grouped_Nff.groups:
                    continue
                Nff_index_to_get = user_grouped_Nff.get_group(Nff).index[0]
                if not dynT_user.ix[Nff_index_to_get].notnull().all():
                    continue
                Nff_index_to_set = Nff_generated_group.index[0]
                if Nff_generated_group['cardinality'][Nff_index_to_set] != highest_cardinality:
                    current_dTT_MS_cohort = most_frequent_axis_dTT_MS_cohort + (Nff_generated_group['id_axis'][Nff_index_to_set] - N_cohort_generated_group['id_axis'][N_cohort_index_to_set]) / (4 * most_frequent_MS_a_cohort) 
                    dynT_tmp_dataframe.ix[Nff_index_to_set]['dTT_MS_cohort'] = current_dTT_MS_cohort                   
                columns_to_set = dynT_user.columns
                dynT_tmp_dataframe.ix[Nff_index_to_set][columns_to_set] = dynT_user.ix[Nff_index_to_get]
                dynT_tmp_dataframe.ix[Nff_index_to_set]['TT_col_break'] = TT_col_break
    else:
        raise Exception('''%s is an unknown user data completeness value. \
The values can be one of %s''', (str(dynT_user_completeness), 
                                 str([DataCompleteness.MIN, 
                                      DataCompleteness.SHORT, 
                                      DataCompleteness.FULL])))
    
    # 4. complete dimT_user
    organ_dim_list = ['L_blade', 'W_blade', 'L_sheath', 'W_sheath', 'L_internode', 'W_internode']
    if dimT_user_completeness == DataCompleteness.MIN:
        index_to_set = dimT_tmp_dataframe[dimT_tmp_dataframe['id_dim']==dimT_tmp_dataframe['id_dim'][0]].index
        for organ_dim in organ_dim_list:
            dimT_tmp_dataframe[organ_dim][index_to_set] = dimT_user[organ_dim][index_to_set]
    elif dimT_user_completeness == DataCompleteness.SHORT:
        user_grouped = dimT_user.groupby('id_axis')
        for generated_id_dim, generated_group in dimT_tmp_dataframe.groupby('id_dim'):
            generated_id_axis = float(str(int(generated_id_dim))[:-2])
            if generated_id_axis not in user_grouped.groups:
                continue
            user_group = user_grouped.get_group(generated_id_axis)
            user_index_phytomer = user_group['index_phytomer'][user_group.index[-1]]
            user_id_dim = float(''.join([str(int(generated_id_axis)), str(int(user_index_phytomer))]))
            if user_id_dim != generated_id_dim:
                continue
            index_to_get = user_group.index
            index_to_set = generated_group.index
            for organ_dim in organ_dim_list:
                dimT_tmp_dataframe[organ_dim][index_to_set] = dimT_user[organ_dim][index_to_get]     
    elif dimT_user_completeness == DataCompleteness.FULL:
        user_grouped = dimT_user.groupby('id_dim')
        for generated_id_dim, generated_group in dimT_tmp_dataframe.groupby('id_dim'):
            if generated_id_dim not in user_grouped.groups:
                continue
            curr_id_dim_user_group_first_index = user_grouped.get_group(generated_id_dim).index[0]
            if not dimT_user.ix[curr_id_dim_user_group_first_index].notnull().all():
                continue
            user_group = user_grouped.get_group(generated_id_dim)
            dimT_tmp_dataframe.ix[generated_group.index] = dimT_user.ix[user_group.index].values
    else:
        raise Exception('''%s is an unknown user data completeness value. \
The values can be one of %s''', (str(dimT_user_completeness), 
                                 str([DataCompleteness.MIN, 
                                      DataCompleteness.SHORT, 
                                      DataCompleteness.FULL]))) 
    
    # 5. second step of the fit process    
    (axeT_dataframe, 
     phenT_abs_dataframe, 
     phenT_dataframe,
     dimT_abs_dataframe, 
     dynT_dataframe, 
     phenT_first_dataframe,
     HS_GL_SSI_T_dataframe, 
     dimT_dataframe) = _gen_adel_input_data_second(axeT_tmp_dataframe, 
                                                    dimT_tmp_dataframe, 
                                                    dynT_tmp_dataframe, 
                                                    GL_number, 
                                                    TT_bolting, 
                                                    TT_flag_leaf_ligulation, 
                                                    delais_TT_stop_del_axis, 
                                                    final_axes_density)
    
    return axeT_dataframe, dimT_dataframe, phenT_dataframe, phenT_abs_dataframe, \
           dimT_abs_dataframe, dynT_dataframe, phenT_first_dataframe, \
           HS_GL_SSI_T_dataframe, tilleringT_dataframe, cohortT_dataframe


def _gen_adel_input_data_first(plant_number, 
                              decide_child_cohort_probabilities, 
                              MS_leaves_number_probabilities,
                              TT_bolting, 
                              TT_flag_leaf_ligulation,
                              final_axes_density,
                              theoretical_cohorts_cardinalities):    
    '''Generate the input data: first step.'''
    # create axeT_tmp
    axeT_tmp_dataframe = axeT.create_axeT_tmp(plant_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities)
    # create dynT_tmp
    dynT_tmp_dataframe = dynT.create_dynT_tmp(axeT_tmp_dataframe['id_phen'].tolist())
    # create dimT_tmp
    dimT_tmp_dataframe = dimT.create_dimT_tmp(dynT_tmp_dataframe)
    # create tilleringT
    tilleringT_dataframe = axeT.create_tilleringT(0, TT_bolting, TT_flag_leaf_ligulation, plant_number, axeT_tmp_dataframe, final_axes_density)
    # create cohortT
    cohortT_dataframe = axeT.create_cohortT(theoretical_cohorts_cardinalities, axeT_tmp_dataframe['id_cohort_axis'])

    return axeT_tmp_dataframe, dimT_tmp_dataframe, dynT_tmp_dataframe, tilleringT_dataframe, cohortT_dataframe


def _gen_adel_input_data_second(axeT_tmp_dataframe, 
                                dimT_user, 
                                dynT_user, 
                                GL_number, 
                                TT_bolting, 
                                TT_flag_leaf_ligulation, 
                                delais_TT_stop_del_axis,
                                final_axes_density):
    '''Generate the input data: second step.'''
    # calculate decimal_elongated_internode_number
    decimal_elongated_internode_number = dynT.calculate_decimal_elongated_internode_number(dimT_user) 
    # create dynT
    dynT_dataframe = dynT.create_dynT(dynT_user, dimT_user, GL_number, decimal_elongated_internode_number)
    # create phenT_abs
    phenT_abs_dataframe = phenT.create_phenT_abs(dynT_dataframe, decimal_elongated_internode_number)
    # create phenT_first
    phenT_first_dataframe = phenT.create_phenT_first(phenT_abs_dataframe)
    # create phenT
    phenT_dataframe = phenT.create_phenT(phenT_abs_dataframe, phenT_first_dataframe)
    # create axeT
    axeT_dataframe = axeT.create_axeT(axeT_tmp_dataframe, phenT_first_dataframe, dynT_dataframe, TT_bolting, TT_flag_leaf_ligulation, delais_TT_stop_del_axis, final_axes_density)
    # create dimT_abs
    dimT_abs_dataframe = dimT.create_dimT_abs(axeT_dataframe, dimT_user, phenT_abs_dataframe, dynT_dataframe)
    # create dimT
    dimT_dataframe = dimT.create_dimT(dimT_abs_dataframe)
    # create HS_GL_SSI_T 
    HS_GL_SSI_T_dataframe = phenT.create_HS_GL_SSI_T(dynT_dataframe)
    
    return axeT_dataframe, phenT_abs_dataframe, phenT_dataframe, \
           dimT_abs_dataframe, dynT_dataframe, phenT_first_dataframe, \
           HS_GL_SSI_T_dataframe, dimT_dataframe
