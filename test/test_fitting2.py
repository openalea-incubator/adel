import random
random.seed('This is an hashable objet to initialize the basic random number generator')
import numpy
import pandas
from openalea.core.path import path
from adel.fit.adel_input_fitting import fit_adel_input_data_first, fit_adel_input_data_second

plant_ids = range(1,101)
cohort_probabilities = {'3': 0.0, '4': 0.900, '5': 0.967, '6': 0.817, '7': 0.083, '8': 0.0, '9': 0.0, '10': 0.0}
main_stem_leaves_number_probability_distribution = {'10': 0.145, '11': 0.818, '12': 0.036, '13': 0.0, '14': 0.0}
secondary_stem_leaves_number_coefficients = {'a_1': 0.9423, 'a_2': 0.555}
emf_1_main_stem_standard_deviation = 30.0
bolting_date = 500 
flowering_date = 1000

expected_fit_adel_input_data_first_dict = {}
linear_expected_fit_adel_input_data_second_dict = {}
bilinear_expected_fit_adel_input_data_second_dict = {}

def setup_module():

    expected_first_axis_table_dataframe = pandas.read_csv(path('data/first_axis_table.csv'))
    expected_first_parameters_dataframe = pandas.read_csv(path('data/first_parameters.csv'))  
    expected_fit_adel_input_data_first_dict.update(axis_table=expected_first_axis_table_dataframe, 
                                                   parameters=expected_first_parameters_dataframe)
    linear_expected_second_axis_table_dataframe = pandas.read_csv(path('data/linear_second_axis_table.csv'))
    linear_expected_absolute_phen_table_dataframe = pandas.read_csv(path('data/linear_second_absolute_phen_table.csv'))
    linear_expected_relative_phen_table_dataframe = pandas.read_csv(path('data/linear_second_relative_phen_table.csv'))
    linear_expected_second_parameters_dataframe = pandas.read_csv(path('data/linear_second_parameters.csv'))
    linear_expected_fit_adel_input_data_second_dict.update(axis_table=linear_expected_second_axis_table_dataframe, 
                                                    absolute_phen_table=linear_expected_absolute_phen_table_dataframe,
                                                    relative_phen_table=linear_expected_relative_phen_table_dataframe,
                                                    parameters=linear_expected_second_parameters_dataframe)
    
    bilinear_expected_second_axis_table_dataframe = pandas.read_csv(path('data/bilinear_second_axis_table.csv'))
    bilinear_expected_absolute_phen_table_dataframe = pandas.read_csv(path('data/bilinear_second_absolute_phen_table.csv'))
    bilinear_expected_relative_phen_table_dataframe = pandas.read_csv(path('data/bilinear_second_relative_phen_table.csv'))
    bilinear_expected_second_parameters_dataframe = pandas.read_csv(path('data/bilinear_second_parameters.csv'))
    bilinear_expected_fit_adel_input_data_second_dict.update(axis_table=bilinear_expected_second_axis_table_dataframe, 
                                                             absolute_phen_table=bilinear_expected_absolute_phen_table_dataframe,
                                                             relative_phen_table=bilinear_expected_relative_phen_table_dataframe,
                                                             parameters=bilinear_expected_second_parameters_dataframe)


def test_fit_adel_input_data_first():
    fit_adel_input_data_first_dict = fit_adel_input_data_first(plant_ids, cohort_probabilities, main_stem_leaves_number_probability_distribution, 
                                                               secondary_stem_leaves_number_coefficients, emf_1_main_stem_standard_deviation, 
                                                               bolting_date, flowering_date)

          
  
    for key in expected_fit_adel_input_data_first_dict.iterkeys():  
#        # uncomment the following lines (until "end") to save the results to csv files. 
#        import tempfile
#        fitting_results_directory = path(tempfile.mkdtemp(suffix='_fitting_results'))
#        test_table_filepath = path(fitting_results_directory/'first_%s.csv' % key)
#        fit_adel_input_data_first_dict[key].to_csv(test_table_filepath, na_rep='NA', index=False)  
#        print 'The results has been saved in %s' % fitting_results_directory
#        # end
        numpy.testing.assert_array_almost_equal(fit_adel_input_data_first_dict[key].values,
                                                expected_fit_adel_input_data_first_dict[key].values)
                
    
def test_fit_adel_input_data_second_linear(): 
    user_parameters_dataframe = pandas.read_csv(path('data/linear_user_parameters.csv'))
    
    fit_adel_input_data_second_dict = fit_adel_input_data_second(expected_fit_adel_input_data_first_dict['axis_table'], 
                                                                 user_parameters_dataframe)
    
    for key in linear_expected_fit_adel_input_data_second_dict.iterkeys():  
#        # uncomment the following lines (until "end") to save the results to csv files. 
#        import tempfile
#        fitting_results_directory = path(tempfile.mkdtemp(suffix='_fitting_results'))
#        test_table_filepath = path(fitting_results_directory/'linear_second_%s.csv' % key)
#        fit_adel_input_data_second_dict[key].to_csv(test_table_filepath, na_rep='NA', index=False)  
#        print 'The results has been saved in %s' % fitting_results_directory
#        # end
        numpy.testing.assert_array_almost_equal(fit_adel_input_data_second_dict[key].values,
                                                linear_expected_fit_adel_input_data_second_dict[key].values)
        
def test_fit_adel_input_data_second_bilinear(): 
    user_parameters_dataframe = pandas.read_csv(path('data/bilinear_user_parameters.csv'))
    
    fit_adel_input_data_second_dict = fit_adel_input_data_second(expected_fit_adel_input_data_first_dict['axis_table'], 
                                                                 user_parameters_dataframe)
    
    for key in bilinear_expected_fit_adel_input_data_second_dict.iterkeys():  
#        # uncomment the following lines (until "end") to save the results to csv files. 
#        import tempfile
#        fitting_results_directory = path(tempfile.mkdtemp(suffix='_fitting_results'))
#        test_table_filepath = path(fitting_results_directory/'bilinear_second_%s.csv' % key)
#        fit_adel_input_data_second_dict[key].to_csv(test_table_filepath, na_rep='NA', index=False)  
#        print 'The results has been saved in %s' % fitting_results_directory
#        # end 
        numpy.testing.assert_array_almost_equal(fit_adel_input_data_second_dict[key].values,
                                                bilinear_expected_fit_adel_input_data_second_dict[key].values)
    

