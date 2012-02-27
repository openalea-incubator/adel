import random
random.seed('This is an hashable objet to initialize the basic random number generator')
import numpy
import pandas
from openalea.core.path import path
from adel.fit.adel_input_fitting import fit_adel_input_data_first, fit_adel_input_data_second

from ConfigParser import ConfigParser

config = ConfigParser()

expected_fit_adel_input_data_first_dict = {}
linear_expected_fit_adel_input_data_second_dict = {}
bilinear_expected_fit_adel_input_data_second_dict = {}

import tempfile
fitting_results_directory = path(tempfile.mkdtemp(suffix='_fitting_results'))

def setup_module():

    config.read(path('data/test_fitting2/default_input_paramters.cfg'))
    
    expected_first_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/first_axis_table.csv'))
    expected_first_parameters_dataframe = pandas.read_csv(path('data/test_fitting2/first_parameters.csv'))  
    expected_fit_adel_input_data_first_dict.update(axis_table=expected_first_axis_table_dataframe, 
                                                   parameters=expected_first_parameters_dataframe)
    linear_expected_second_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_axis_table.csv'))
    linear_expected_absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_absolute_phen_table.csv'))
    linear_expected_relative_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_relative_phen_table.csv'))
    linear_expected_second_parameters_dataframe = pandas.read_csv(path('data/test_fitting2/linear_second_parameters.csv'))
    linear_expected_fit_adel_input_data_second_dict.update(axis_table=linear_expected_second_axis_table_dataframe, 
                                                    absolute_phen_table=linear_expected_absolute_phen_table_dataframe,
                                                    relative_phen_table=linear_expected_relative_phen_table_dataframe,
                                                    parameters=linear_expected_second_parameters_dataframe)
    
    bilinear_expected_second_axis_table_dataframe = pandas.read_csv(path('data/test_fitting2/bilinear_second_axis_table.csv'))
    bilinear_expected_absolute_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/bilinear_second_absolute_phen_table.csv'))
    bilinear_expected_relative_phen_table_dataframe = pandas.read_csv(path('data/test_fitting2/bilinear_second_relative_phen_table.csv'))
    bilinear_expected_second_parameters_dataframe = pandas.read_csv(path('data/test_fitting2/bilinear_second_parameters.csv'))
    bilinear_expected_fit_adel_input_data_second_dict.update(axis_table=bilinear_expected_second_axis_table_dataframe, 
                                                             absolute_phen_table=bilinear_expected_absolute_phen_table_dataframe,
                                                             relative_phen_table=bilinear_expected_relative_phen_table_dataframe,
                                                             parameters=bilinear_expected_second_parameters_dataframe)                                              


def test_fit_adel_input_data_first():

    cohort_probabilities = dict([(name_value_tuple[0], float(name_value_tuple[1])) for name_value_tuple in config.items('cohort_probabilities')])
    main_stem_leaves_number_probability_distribution = dict([(name_value_tuple[0], float(name_value_tuple[1])) for name_value_tuple in config.items('main_stem_leaves_number_probability_distribution')])

    fit_adel_input_data_first_dict = fit_adel_input_data_first(plant_number=int(config._sections['plant_number']['plant_number']), 
                                                               cohort_probabilities=cohort_probabilities, 
                                                               main_stem_leaves_number_probability_distribution=main_stem_leaves_number_probability_distribution,      
                                                               bolting_date=int(config._sections['bolting_date']['bolting_date']), 
                                                               flowering_date=int(config._sections['flowering_date']['flowering_date']))
          
    for key in expected_fit_adel_input_data_first_dict.iterkeys():  
        test_table_filepath = path(fitting_results_directory/'first_%s.csv' % key)
        fit_adel_input_data_first_dict[key].to_csv(test_table_filepath, na_rep='NA', index=False)  
        print 'The results has been saved in %s' % fitting_results_directory
        #uncomment the 2 next lines to open the csv result file with the default editor. ATTENTION: for windows only.
        #import os
        #os.system("start %s" % test_table_filepath)
        numpy.testing.assert_array_almost_equal(fit_adel_input_data_first_dict[key].values,
                                                expected_fit_adel_input_data_first_dict[key].values)
                
    
def test_fit_adel_input_data_second_linear(): 
    user_parameters_dataframe = pandas.read_csv(path('data/test_fitting2/linear_user_parameters.csv'))
    
    fit_adel_input_data_second_dict = fit_adel_input_data_second(expected_fit_adel_input_data_first_dict['axis_table'], 
                                                                 user_parameters_dataframe)
    
    for key in linear_expected_fit_adel_input_data_second_dict.iterkeys():  
        test_table_filepath = path(fitting_results_directory/'linear_second_%s.csv' % key)
        fit_adel_input_data_second_dict[key].to_csv(test_table_filepath, na_rep='NA', index=False)  
        print 'The results has been saved in %s' % fitting_results_directory
        #uncomment the 2 next lines to open the csv result file with the default editor. ATTENTION: for windows only.
        #import os
        #os.system("start %s" % test_table_filepath)
        numpy.testing.assert_array_almost_equal(fit_adel_input_data_second_dict[key].values,
                                                linear_expected_fit_adel_input_data_second_dict[key].values)
        
def test_fit_adel_input_data_second_bilinear(): 
    user_parameters_dataframe = pandas.read_csv(path('data/test_fitting2/bilinear_user_parameters.csv'))
    
    fit_adel_input_data_second_dict = fit_adel_input_data_second(expected_fit_adel_input_data_first_dict['axis_table'], 
                                                                 user_parameters_dataframe)
    
    for key in bilinear_expected_fit_adel_input_data_second_dict.iterkeys():  
        test_table_filepath = path(fitting_results_directory/'bilinear_second_%s.csv' % key)
        fit_adel_input_data_second_dict[key].to_csv(test_table_filepath, na_rep='NA', index=False)  
        print 'The results has been saved in %s' % fitting_results_directory
        #uncomment the 2 next lines to open the csv result file with the default editor. ATTENTION: for windows only.
        #import os
        #os.system("start %s" % test_table_filepath)
        numpy.testing.assert_array_almost_equal(fit_adel_input_data_second_dict[key].values,
                                                bilinear_expected_fit_adel_input_data_second_dict[key].values)
    

