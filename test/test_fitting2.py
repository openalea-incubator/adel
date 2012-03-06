import random
from adel.fit import axis_table_fitting, dim_table_fitting, phen_table_fitting, fit_adel_input_data_first, fit_adel_input_data_second
random.seed('This is an hashable objet to initialize the basic random number generator')
import numpy
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

import tempfile
fitting_results_directory = path(tempfile.mkdtemp(suffix='_fitting_results'))

        
def test_fit_axis_table_first():
    expected_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/random_first_axis_table.csv'))
    axis_table_dataframe = axis_table_fitting.fit_axis_table_first(plant_number, cohort_probabilities, main_stem_leaves_number_probability_distribution, bolting_date, flowering_date)
    test_table_filepath = path(fitting_results_directory/'linear_first_axis_table.csv')
    axis_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    numpy.testing.assert_array_almost_equal(axis_table_dataframe.values,
                                            expected_axis_table_dataframe.values)


def test_fit_user_parameters_first():
    axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_axis_table.csv'))
    expected_parameters_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_parameters_table.csv'))
    parameters_table_dataframe = fit_adel_input_data_first.fit_user_parameters_first(axis_table_dataframe['id_phen'].tolist())
    test_table_filepath = path(fitting_results_directory/'linear_first_parameters_table.csv')
    parameters_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    numpy.testing.assert_array_almost_equal(parameters_table_dataframe.values,
                                            expected_parameters_table_dataframe.values)


def test_fit_dim_table_first():
    parameters_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_parameters_table.csv'))
    expected_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_dim_table.csv'))
    dim_table_dataframe = dim_table_fitting.fit_dim_table_first(parameters_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'linear_first_dim_table.csv')
    dim_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    numpy.testing.assert_array_almost_equal(dim_table_dataframe.values,
                                            expected_dim_table_dataframe.values)

        
def test_fit_user_parameters_second_linear():
    user_parameter_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_user_parameters_table.csv'))
    expected_fitted_parameter_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_parameters_table.csv'))
    second_parameters_dataframe = fit_adel_input_data_second.fit_user_parameters_second(user_parameter_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'linear_second_parameters_table.csv')
    second_parameters_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    numpy.testing.assert_array_almost_equal(second_parameters_dataframe.values,
                                            expected_fitted_parameter_dataframe.values)
    

def test_fit_phen_table_second_linear():
    second_parameters_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_parameters_table.csv'))
    expected_absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_absolute_phen_table.csv'))
    absolute_phen_table_dataframe = phen_table_fitting.fit_phen_table_second(second_parameters_dataframe)
    test_table_filepath = path(fitting_results_directory/'linear_second_absolute_phen_table.csv')
    absolute_phen_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    numpy.testing.assert_array_almost_equal(absolute_phen_table_dataframe.values,
                                            expected_absolute_phen_table_dataframe.values)


def test_create_phen_table_relative_dataframe_linear():
    absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_absolute_phen_table.csv'))
    expected_relative_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_relative_phen_table.csv'))
    relative_phen_table_dataframe = phen_table_fitting.create_phen_table_relative_dataframe(absolute_phen_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'linear_second_relative_phen_table.csv')
    relative_phen_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    numpy.testing.assert_array_almost_equal(relative_phen_table_dataframe.values,
                                            expected_relative_phen_table_dataframe.values)


def test_fit_axis_table_second_linear():
    first_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_first_axis_table.csv'))
    expected_second_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_axis_table.csv'))
    second_axis_table_dataframe = axis_table_fitting.fit_axis_table_second(first_axis_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'linear_second_axis_table.csv')
    second_axis_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    numpy.testing.assert_array_almost_equal(second_axis_table_dataframe.values,
                                            expected_second_axis_table_dataframe.values)
 

def test_fit_dim_table_second_linear():
    user_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_user_dim_table.csv'))
    absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_absolute_phen_table.csv'))
    expected_dim_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_dim_table.csv'))
    dim_table_dataframe = dim_table_fitting.fit_dim_table_second(user_dim_table_dataframe, absolute_phen_table_dataframe)
    test_table_filepath = path(fitting_results_directory/'linear_second_dim_table.csv')
    dim_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    # ATTENTION: for windows only.
    # Uncomment the 2 next lines to open the csv result file with the default editor. 
    #import os
    #os.system("start %s" % test_table_filepath)
    numpy.testing.assert_array_almost_equal(dim_table_dataframe.values,
                                            expected_dim_table_dataframe.values)

