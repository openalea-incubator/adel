import random
from adel.fit import axis_table_fitting, dim_table_fitting, phen_table_fitting, fit_adel_input_data_first, fit_adel_input_data_second
random.seed(1234)
import numpy as np
import pandas
from openalea.core.path import path

from ConfigParser import ConfigParser

config = ConfigParser()
config.read(path('data/test_fitting2/default_input_parameters.cfg'))
plant_number = int(config._sections['plant_number']['plant_number'])
cohort_probabilities = dict([(name_value_tuple[0], float(name_value_tuple[1])) for name_value_tuple in config.items('cohort_probabilities')])
main_stem_leaves_number_probability_distribution = dict([(name_value_tuple[0], float(name_value_tuple[1])) for name_value_tuple in config.items('main_stem_leaves_number_probability_distribution')])
bolting_date = int(config._sections['bolting_date']['bolting_date'])
flowering_date = int(config._sections['flowering_date']['flowering_date'])
GL_number = dict([(float(name_value_tuple[0]), float(name_value_tuple[1])) for name_value_tuple in config.items('GL_number')])
delais_TT_stop_del_axis = int(config._sections['delais_TT_stop_del_axis']['delais_tt_stop_del_axis'])
final_axes_number = int(config._sections['final_axes_number']['final_axes_number'])

import tempfile
fitting_results_directory = path(tempfile.mkdtemp(suffix='_fitting_results'))
relative_tolerance = 10e-3
absolute_tolerance = 10e-3
#TODO: check with Mariem that is sufficient
        
def test_fit_axis_table_first():
    expected_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/random_first_axis_table.csv'))
    axis_table_dataframe = axis_table_fitting.fit_axis_table_first(plant_number, cohort_probabilities, main_stem_leaves_number_probability_distribution)
    test_table_filepath = path(fitting_results_directory/'result_linear_first_axis_table.csv')
    axis_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(axis_table_dataframe.values, expected_axis_table_dataframe.values, relative_tolerance, absolute_tolerance)


def test_fit_user_parameters_first():
    axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_axis_table.csv'))
    expected_parameters_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_parameters_table.csv'))
    parameters_table_dataframe = fit_adel_input_data_first.fit_user_parameters_first(axis_table_dataframe['id_phen'].tolist())
    test_table_filepath = path(fitting_results_directory/'result_linear_first_parameters_table.csv')
    parameters_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(parameters_table_dataframe.values, expected_parameters_table_dataframe.values, relative_tolerance, absolute_tolerance)


def test_fit_dim_table_first():
    parameters_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_parameters_table.csv'))
    expected_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_dim_table.csv'))
    dim_table_dataframe = dim_table_fitting.fit_dim_table_first(parameters_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'result_linear_first_dim_table.csv')
    dim_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(dim_table_dataframe.values, expected_dim_table_dataframe.values, relative_tolerance, absolute_tolerance)

        
def test_fit_user_parameters_second_linear():
    user_parameter_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_user_parameters_table.csv'))
    user_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_user_dim_table.csv'))
    expected_fitted_parameter_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_parameters_table.csv'))
    second_parameters_dataframe = fit_adel_input_data_second.fit_user_parameters_second(user_parameter_table_dataframe, user_dim_table_dataframe, GL_number)
    test_table_filepath = path(fitting_results_directory/'result_linear_second_parameters_table.csv')
    second_parameters_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(second_parameters_dataframe.values, expected_fitted_parameter_dataframe.values, relative_tolerance, absolute_tolerance)
    

def test_fit_phen_table_second_linear():
    second_parameters_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_parameters_table.csv'))
    expected_absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_absolute_phen_table.csv'))
    absolute_phen_table_dataframe = phen_table_fitting.fit_phen_table_second(second_parameters_dataframe)
    test_table_filepath = path(fitting_results_directory/'result_linear_second_absolute_phen_table.csv')
    absolute_phen_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(absolute_phen_table_dataframe.values, expected_absolute_phen_table_dataframe.values, relative_tolerance, absolute_tolerance)


def test_create_first_leaf_phen_table_dataframe_linear():
    absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_absolute_phen_table.csv'))
    expected_first_leaf_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_leaf_phen_table.csv'))
    first_leaf_phen_table_dataframe = phen_table_fitting.create_first_leaf_phen_table_dataframe(absolute_phen_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'result_first_leaf_phen_table.csv')
    first_leaf_phen_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(first_leaf_phen_table_dataframe.values, expected_first_leaf_phen_table_dataframe.values, relative_tolerance, absolute_tolerance)


def test_create_phen_table_relative_dataframe_linear():
    absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_absolute_phen_table.csv'))
    expected_relative_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_relative_phen_table.csv'))
    first_leaf_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_leaf_phen_table.csv'))
    relative_phen_table_dataframe = phen_table_fitting.create_phen_table_relative_dataframe(absolute_phen_table_dataframe, first_leaf_phen_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'result_linear_second_relative_phen_table.csv')
    relative_phen_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(relative_phen_table_dataframe.values, expected_relative_phen_table_dataframe.values, relative_tolerance, absolute_tolerance)


def test_create_HS_GL_SSI_dynamic_dataframe_linear():
    second_parameters_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_parameters_table.csv'))
    expected_HS_GL_SSI_dynamic_dataframe = pandas.read_csv(path('data/test_fitting2/linear_HS_GL_SSI_dynamic_table.csv'))
    HS_GL_SSI_dynamic_dataframe = phen_table_fitting.create_HS_GL_SSI_dynamic_dataframe(second_parameters_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'result_HS_GL_SSI_dynamic_table.csv')
    HS_GL_SSI_dynamic_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(HS_GL_SSI_dynamic_dataframe.values, expected_HS_GL_SSI_dynamic_dataframe.values, relative_tolerance, absolute_tolerance)


def test_fit_axis_table_second_linear():
    first_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_axis_table.csv'))
    expected_second_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_axis_table.csv'))
    first_leaf_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_leaf_phen_table.csv'))
    second_axis_table_dataframe = axis_table_fitting.fit_axis_table_second(first_axis_table_dataframe, first_leaf_phen_table_dataframe, bolting_date, flowering_date, delais_TT_stop_del_axis, final_axes_number)
    test_table_filepath = path(fitting_results_directory/'result_linear_second_axis_table.csv')
    second_axis_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(second_axis_table_dataframe.values, expected_second_axis_table_dataframe.values, relative_tolerance, absolute_tolerance)
 

def test_fit_dim_table_second_linear():
    user_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_user_dim_table.csv'))
    absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_absolute_phen_table.csv'))
    expected_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_dim_table.csv'))
    dim_table_dataframe = dim_table_fitting.fit_dim_table_second(user_dim_table_dataframe, absolute_phen_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'result_linear_second_dim_table.csv')
    dim_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(dim_table_dataframe.values, expected_dim_table_dataframe.values, relative_tolerance, absolute_tolerance)


def test_create_tillering_dynamic_dataframe():
    axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_axis_table.csv'))
    expected_tillering_dynamic_dataframe = pandas.read_csv(path('data/test_fitting2/tillering_dynamic_table.csv'))
    tillering_dynamic_dataframe = axis_table_fitting.create_tillering_dynamic_dataframe(0, bolting_date, flowering_date, plant_number, axis_table_dataframe, final_axes_number)
    test_table_filepath = path(fitting_results_directory/'result_tillering_dynamic_table.csv')
    tillering_dynamic_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(tillering_dynamic_dataframe.values, expected_tillering_dynamic_dataframe.values, relative_tolerance, absolute_tolerance)


def test_create_dim_table_relative_dataframe_linear():
    absolute_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_dim_table.csv'))
    expected_relative_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_relative_dim_table.csv'))
    relative_dim_table_dataframe = dim_table_fitting.create_dim_table_relative_dataframe(absolute_dim_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'result_linear_relative_dim_table.csv')
    relative_dim_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    np.testing.assert_allclose(relative_dim_table_dataframe.values, expected_relative_dim_table_dataframe.values, relative_tolerance, absolute_tolerance)


if __name__ == '__main__':
    test_fit_axis_table_first()
    test_fit_user_parameters_first()
    test_fit_dim_table_first()
    test_fit_user_parameters_second_linear()
    test_fit_phen_table_second_linear()
    test_create_first_leaf_phen_table_dataframe_linear()
    test_create_phen_table_relative_dataframe_linear()
    test_create_HS_GL_SSI_dynamic_dataframe_linear()
    test_fit_axis_table_second_linear()
    test_create_tillering_dynamic_dataframe()
    test_fit_dim_table_second_linear()
    test_create_dim_table_relative_dataframe_linear()
