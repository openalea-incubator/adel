# -*- python -*-
#
#       Adel.Fit
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
'''This module is a facade for the second step of Adel input data fitting .
'''  
from adel.fit import axis_table_fitting, phen_table_fitting, organ_dimensions_table_fitting, leaf_dynamic_parameters_table_fitting

def fit_adel_input_data_second(first_axis_table_dataframe, user_organ_dimensions_table_dataframe, user_parameter_table_dataframe, 
                               GL_number={1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}, 
                               bolting_date=500, flowering_date=1440, 
                               delais_TT_stop_del_axis=600,
                               final_axes_number=250):
    '''
    Complete the axis table, the dim table and the leaf_dynamic_parameters table, and create the absolute/relative phen tables.
    Construct:
        - the parameter table which contains the information about tillers behaviour. The leaf_dynamic_parameters table is not an input of ADEL.
          It is used only to build the other tables: ParametersTable
        - ADEL input data tables: second_axis_table_dataframe, relative_second_phen_table_dataframe, relative_organ_dimensions_table_dataframe
        - tables in order the user could check the fitting results: first_leaf_phen_table_dataframe, 
          HS_GL_SSI_dynamic_dataframe
    
    :Parameters:
        - `first_axis_table_dataframe` : the first axis table.
        - `user_organ_dimensions_table_dataframe` : the user dim table.
        - `user_parameter_table_dataframe` : the user leaf_dynamic_parameters table.
        - `GL_number` : the GL decimal number measured at several thermal time (including the senescence end).
        - `bolting_date` : The bolting date. Must be positive or null, and lesser than flowering_date.
        - `flowering_date` : The flowering date. Must be positive or null, and greater than bolting_date.
        - `delais_TT_stop_del_axis` : Thermal time during which a tiller remains present on the plant after the tiller has stopped 
                                      growing (for tillers that senesce).
        - `final_axes_number` : the final number of axes per square meter.
    :Types:
        - `first_axis_table_dataframe` : pandas.DataFrame
        - `user_organ_dimensions_table_dataframe` : pandas.DataFrame
        - `user_parameter_table_dataframe` : pandas.DataFrame
        - `GL_number` : a dict of 2-tuples
        - `bolting_date` : int
        - `flowering_date` : int
        - `delais_TT_stop_del_axis` : int
        - `final_axes_number` : int

    :return: The fitted axis table, the fitted absolute phen table, the fitted relative phen table,
    the fitted dim table and the fitted leaf_dynamic_parameters table.
    :rtype: a tuple of pandas.DataFrame
    '''
    # Fit the leaf_dynamic_parameters provided by the user
    second_leaf_dynamic_parameters_table_dataframe = leaf_dynamic_parameters_table_fitting.fit_user_leaf_dynamic_parameters_second(user_parameter_table_dataframe, user_organ_dimensions_table_dataframe, GL_number)
    # Fit absolute PhenTable
    absolute_second_phen_table_dataframe = phen_table_fitting.fit_phen_table_second(second_leaf_dynamic_parameters_table_dataframe)
    # Extract the first leaves data from absolute_second_phen_table_dataframe
    first_leaf_phen_table_dataframe = phen_table_fitting.create_first_leaf_phen_table_dataframe(absolute_second_phen_table_dataframe)
    # Fit relative PhenTable
    relative_second_phen_table_dataframe = phen_table_fitting.create_phen_table_relative_dataframe(absolute_second_phen_table_dataframe, first_leaf_phen_table_dataframe)
    # Fit AxisTable
    second_axis_table_dataframe = axis_table_fitting.fit_axis_table_second(first_axis_table_dataframe, first_leaf_phen_table_dataframe, bolting_date, flowering_date, delais_TT_stop_del_axis, final_axes_number)
    # Fit DimTable
    absolute_organ_dimensions_table_dataframe = organ_dimensions_table_fitting.fit_organ_dimensions_table_second(user_organ_dimensions_table_dataframe, absolute_second_phen_table_dataframe)
    # Fit relative dimTable
    relative_organ_dimensions_table_dataframe = organ_dimensions_table_fitting.create_organ_dimensions_table_relative_dataframe(absolute_organ_dimensions_table_dataframe)
    # Create a table with the following columns: id_axis,TT,HS,GL,SSI
    HS_GL_SSI_dynamic_dataframe = phen_table_fitting.create_HS_GL_SSI_dynamic_dataframe(second_leaf_dynamic_parameters_table_dataframe)
    
    return second_axis_table_dataframe, absolute_second_phen_table_dataframe, relative_second_phen_table_dataframe, \
           absolute_organ_dimensions_table_dataframe, second_leaf_dynamic_parameters_table_dataframe, first_leaf_phen_table_dataframe, \
           HS_GL_SSI_dynamic_dataframe, relative_organ_dimensions_table_dataframe
    
    
    
    
