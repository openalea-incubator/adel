import random

from adel.fit import axis_table_fitting, dim_table_fitting, phen_table_fitting, fit_adel_input_data_first, fit_adel_input_data_second,\
    fit_adel_input_data, parameters_table_fitting
import numpy as np
import pandas
from openalea.core.path import path
from openalea.core.alea import *
from nose.tools import with_setup

random.seed(1234)

initial_random_state = random.getstate()

plant_number = 100
cohort_probabilities = {'3': 0.0, '4': 0.900, '5': 0.967, '6': 0.817, '7': 0.083, '8': 0.0, '9': 0.0, '10': 0.0}
main_stem_leaves_number_probability_distribution = {'10': 0.145, '11': 0.818, '12': 0.036, '13': 0.0, '14': 0.0}
bolting_date = 500
flowering_date = 1440
final_axes_number = 250
GL_number = {1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}
delais_TT_stop_del_axis = 600
dTT_MS_cohort = {'4': 70, '5': 80, '6': 90, '7': 100, '8': 110, '9': 120, '10': 130}

expected_results_dir = path('data/test_fitting2')
default_expected_results_dir = expected_results_dir.joinpath('default')
min_min_expected_results_dir = expected_results_dir.joinpath('min_min')
short_short_expected_results_dir = expected_results_dir.joinpath('short_short')
full_full_expected_results_dir = expected_results_dir.joinpath('full_full')

import tempfile
fitting_results_directory = path(tempfile.mkdtemp(suffix='_fitting_results'))
default_results = fitting_results_directory.joinpath('default')
if not default_results.exists():
    default_results.mkdir()
relative_tolerance = 10e-3
absolute_tolerance = 10e-3


def reinit_random_state():
    random.setstate(initial_random_state)


@with_setup(reinit_random_state)
def test_generate_axes():
    expected_axis_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_axis_table.csv')
    axis_table_dataframe = axis_table_fitting.generate_axes(plant_number, cohort_probabilities, main_stem_leaves_number_probability_distribution)
    test_table_filepath = default_results.joinpath('linear_first_axis_table.csv')
    axis_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(axis_table_dataframe.values, expected_axis_table_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_fit_user_parameters_first():
    axis_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_axis_table.csv')
    expected_parameters_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_parameters_table.csv')
    parameters_table_dataframe = parameters_table_fitting.fit_user_parameters_first(axis_table_dataframe['id_phen'].tolist())
    test_table_filepath = default_results.joinpath('linear_first_parameters_table.csv')
    parameters_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(parameters_table_dataframe.values, expected_parameters_table_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_fit_dim_table_first():
    parameters_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_parameters_table.csv')
    expected_dim_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_dim_table.csv')
    dim_table_dataframe = dim_table_fitting.fit_dim_table_first(parameters_table_dataframe)
    test_table_filepath = default_results.joinpath('linear_first_dim_table.csv')
    dim_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dim_table_dataframe.values, expected_dim_table_dataframe.values, relative_tolerance, absolute_tolerance)

        
@with_setup(reinit_random_state)
def test_fit_user_parameters_second_linear():
    user_parameter_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_user_parameters_table.csv')
    user_dim_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_user_dim_table.csv')
    expected_fitted_parameter_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_parameters_table.csv')
    second_parameters_dataframe = parameters_table_fitting.fit_user_parameters_second(user_parameter_table_dataframe, user_dim_table_dataframe, GL_number)
    test_table_filepath = default_results.joinpath('linear_second_parameters_table.csv')
    second_parameters_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(second_parameters_dataframe.values, expected_fitted_parameter_dataframe.values, relative_tolerance, absolute_tolerance)
    

@with_setup(reinit_random_state)
def test_fit_phen_table_second_linear():
    second_parameters_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_parameters_table.csv')
    expected_absolute_phen_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_absolute_phen_table.csv')
    absolute_phen_table_dataframe = phen_table_fitting.fit_phen_table_second(second_parameters_dataframe)
    test_table_filepath = default_results.joinpath('linear_second_absolute_phen_table.csv')
    absolute_phen_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(absolute_phen_table_dataframe.values, expected_absolute_phen_table_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_first_leaf_phen_table_dataframe_linear():
    absolute_phen_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_absolute_phen_table.csv')
    expected_first_leaf_phen_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_leaf_phen_table.csv')
    first_leaf_phen_table_dataframe = phen_table_fitting.create_first_leaf_phen_table_dataframe(absolute_phen_table_dataframe)
    test_table_filepath = default_results.joinpath('linear_first_leaf_phen_table.csv')
    first_leaf_phen_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(first_leaf_phen_table_dataframe.values, expected_first_leaf_phen_table_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_phen_table_relative_dataframe_linear():
    absolute_phen_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_absolute_phen_table.csv')
    expected_relative_phen_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_relative_phen_table.csv')
    first_leaf_phen_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_leaf_phen_table.csv')
    relative_phen_table_dataframe = phen_table_fitting.create_phen_table_relative_dataframe(absolute_phen_table_dataframe, first_leaf_phen_table_dataframe)
    test_table_filepath = default_results.joinpath('linear_second_relative_phen_table.csv')
    relative_phen_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(relative_phen_table_dataframe.values, expected_relative_phen_table_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_HS_GL_SSI_dynamic_dataframe_linear():
    second_parameters_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_parameters_table.csv')
    expected_HS_GL_SSI_dynamic_dataframe = pandas.read_csv(default_expected_results_dir/'linear_HS_GL_SSI_dynamic_table.csv')
    HS_GL_SSI_dynamic_dataframe = phen_table_fitting.create_HS_GL_SSI_dynamic_dataframe(second_parameters_table_dataframe)
    test_table_filepath = default_results.joinpath('linear_HS_GL_SSI_dynamic_table.csv')
    HS_GL_SSI_dynamic_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(HS_GL_SSI_dynamic_dataframe.values, expected_HS_GL_SSI_dynamic_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_fit_axis_table_second_linear():
    first_axis_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_axis_table.csv')
    expected_second_axis_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_axis_table.csv')
    first_leaf_phen_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_leaf_phen_table.csv')
    second_axis_table_dataframe = axis_table_fitting.fit_axis_table_second(first_axis_table_dataframe, first_leaf_phen_table_dataframe, bolting_date, flowering_date, delais_TT_stop_del_axis, final_axes_number)
    test_table_filepath = default_results.joinpath('linear_second_axis_table.csv')
    second_axis_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(second_axis_table_dataframe.values, expected_second_axis_table_dataframe.values, relative_tolerance, absolute_tolerance)
 

@with_setup(reinit_random_state)
def test_fit_dim_table_second_linear():
    user_dim_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_user_dim_table.csv')
    absolute_phen_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_absolute_phen_table.csv')
    expected_dim_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_dim_table.csv')
    dim_table_dataframe = dim_table_fitting.fit_dim_table_second(user_dim_table_dataframe, absolute_phen_table_dataframe)
    test_table_filepath = default_results.joinpath('linear_second_dim_table.csv')
    dim_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dim_table_dataframe.values, expected_dim_table_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_tillering_dynamic_dataframe():
    axis_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_first_axis_table.csv')
    expected_tillering_dynamic_dataframe = pandas.read_csv(default_expected_results_dir/'tillering_dynamic_table.csv')
    tillering_dynamic_dataframe = axis_table_fitting.create_tillering_dynamic_dataframe(0, bolting_date, flowering_date, plant_number, axis_table_dataframe, final_axes_number)
    test_table_filepath = default_results.joinpath('tillering_dynamic_table.csv')
    tillering_dynamic_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(tillering_dynamic_dataframe.values, expected_tillering_dynamic_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_dim_table_relative_dataframe_linear():
    absolute_dim_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_second_dim_table.csv')
    expected_relative_dim_table_dataframe = pandas.read_csv(default_expected_results_dir/'linear_relative_dim_table.csv')
    relative_dim_table_dataframe = dim_table_fitting.create_dim_table_relative_dataframe(absolute_dim_table_dataframe)
    test_table_filepath = default_results.joinpath('linear_relative_dim_table.csv')
    relative_dim_table_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(relative_dim_table_dataframe.values, expected_relative_dim_table_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_fit_adel_input_data_min_min():
    user_parameters_completeness = fit_adel_input_data.DataCompleteness.MIN
    user_dims_completeness = fit_adel_input_data.DataCompleteness.MIN
    TT_col_break = 0.0
    user_parameters = {'a_cohort': 0.0102, 
                       'TT_col_0': -0.771289027, 
                       'TT_col_nff': 1078.0, 
                       'n0': 4.871559739,
                       'n1': 3.24283148,
                       'n2': 5.8}
    user_dims = pandas.read_csv(min_min_expected_results_dir/'user_dim_table.csv')
    results = fit_adel_input_data.fit_adel_input_data(user_parameters,
                                                      user_dims, 
                                                      plant_number, 
                                                      cohort_probabilities, 
                                                      main_stem_leaves_number_probability_distribution, 
                                                      bolting_date, 
                                                      flowering_date, 
                                                      final_axes_number, 
                                                      GL_number, 
                                                      delais_TT_stop_del_axis, 
                                                      dTT_MS_cohort, 
                                                      TT_col_break,
                                                      user_parameters_completeness, 
                                                      user_dims_completeness)
    
    expected_axis_table = pandas.read_csv(min_min_expected_results_dir/'linear_second_axis_table.csv')
    expected_absolute_phen_table = pandas.read_csv(min_min_expected_results_dir/'linear_second_absolute_phen_table.csv')
    expected_relative_phen_table = pandas.read_csv(min_min_expected_results_dir/'linear_second_relative_phen_table.csv')
    expected_absolute_dim_table = pandas.read_csv(min_min_expected_results_dir/'linear_second_dim_table.csv')
    expected_parameters_table = pandas.read_csv(min_min_expected_results_dir/'linear_second_parameters_table.csv')
    expected_first_leaf_phen_table = pandas.read_csv(min_min_expected_results_dir/'linear_first_leaf_phen_table.csv')
    expected_HS_GL_SSI_dynamic_table = pandas.read_csv(min_min_expected_results_dir/'linear_HS_GL_SSI_dynamic_table.csv')
    expected_relative_dim_table = pandas.read_csv(min_min_expected_results_dir/'linear_relative_dim_table.csv')
    expected_tillering_dynamic_table = pandas.read_csv(min_min_expected_results_dir/'tillering_dynamic_table.csv')

    to_compare = {'linear_second_axis_table': (expected_axis_table, results[0]),
                 'linear_second_absolute_phen_table': (expected_absolute_phen_table, results[1]),
                 'linear_second_relative_phen_table': (expected_relative_phen_table, results[2]),
                 'linear_second_dim_table': (expected_absolute_dim_table, results[3]),
                 'linear_second_parameters_table': (expected_parameters_table, results[4]),
                 'linear_first_leaf_phen_table': (expected_first_leaf_phen_table, results[5]),
                 'linear_HS_GL_SSI_dynamic_table': (expected_HS_GL_SSI_dynamic_table, results[6]),
                 'linear_relative_dim_table': (expected_relative_dim_table, results[7]),
                 'tillering_dynamic_table': (expected_tillering_dynamic_table, results[8])}
    
    _check_results(to_compare, user_parameters_completeness, user_dims_completeness)

         
@with_setup(reinit_random_state)
def test_fit_adel_input_data_short_short():
    user_parameters_completeness = fit_adel_input_data.DataCompleteness.SHORT
    user_dims_completeness = fit_adel_input_data.DataCompleteness.SHORT
    TT_col_break = 0.0
    user_parameters = pandas.read_csv(short_short_expected_results_dir/'user_parameters_table.csv')
    user_dims = pandas.read_csv(short_short_expected_results_dir/'user_dim_table.csv')
    results = fit_adel_input_data.fit_adel_input_data(user_parameters,
                                                      user_dims, 
                                                      plant_number, 
                                                      cohort_probabilities, 
                                                      main_stem_leaves_number_probability_distribution, 
                                                      bolting_date, 
                                                      flowering_date, 
                                                      final_axes_number, 
                                                      GL_number, 
                                                      delais_TT_stop_del_axis, 
                                                      dTT_MS_cohort, 
                                                      TT_col_break, 
                                                      user_parameters_completeness, 
                                                      user_dims_completeness)
    
    expected_axis_table = pandas.read_csv(short_short_expected_results_dir/'linear_second_axis_table.csv')
    expected_absolute_phen_table = pandas.read_csv(short_short_expected_results_dir/'linear_second_absolute_phen_table.csv')
    expected_relative_phen_table = pandas.read_csv(short_short_expected_results_dir/'linear_second_relative_phen_table.csv')
    expected_absolute_dim_table = pandas.read_csv(short_short_expected_results_dir/'linear_second_dim_table.csv')
    expected_parameters_table = pandas.read_csv(short_short_expected_results_dir/'linear_second_parameters_table.csv')
    expected_first_leaf_phen_table = pandas.read_csv(short_short_expected_results_dir/'linear_first_leaf_phen_table.csv')
    expected_HS_GL_SSI_dynamic_table = pandas.read_csv(short_short_expected_results_dir/'linear_HS_GL_SSI_dynamic_table.csv')
    expected_relative_dim_table = pandas.read_csv(short_short_expected_results_dir/'linear_relative_dim_table.csv')
    expected_tillering_dynamic_table = pandas.read_csv(short_short_expected_results_dir/'tillering_dynamic_table.csv')

    to_compare = {'linear_second_axis_table': (expected_axis_table, results[0]),
                 'linear_second_absolute_phen_table': (expected_absolute_phen_table, results[1]),
                 'linear_second_relative_phen_table': (expected_relative_phen_table, results[2]),
                 'linear_second_dim_table': (expected_absolute_dim_table, results[3]),
                 'linear_second_parameters_table': (expected_parameters_table, results[4]),
                 'linear_first_leaf_phen_table': (expected_first_leaf_phen_table, results[5]),
                 'linear_HS_GL_SSI_dynamic_table': (expected_HS_GL_SSI_dynamic_table, results[6]),
                 'linear_relative_dim_table': (expected_relative_dim_table, results[7]),
                 'tillering_dynamic_table': (expected_tillering_dynamic_table, results[8])}
    
    _check_results(to_compare, user_parameters_completeness, user_dims_completeness)


@with_setup(reinit_random_state)
def test_fit_adel_input_data_full_full():
    user_parameters_completeness = fit_adel_input_data.DataCompleteness.FULL
    user_dims_completeness = fit_adel_input_data.DataCompleteness.FULL
    TT_col_break = 0.0
    user_parameters = pandas.read_csv(full_full_expected_results_dir/'user_parameters_table.csv')
    user_dims = pandas.read_csv(full_full_expected_results_dir/'user_dim_table.csv')
    results = fit_adel_input_data.fit_adel_input_data(user_parameters,
                                                      user_dims,
                                                      plant_number, 
                                                      cohort_probabilities, 
                                                      main_stem_leaves_number_probability_distribution, 
                                                      bolting_date, 
                                                      flowering_date, 
                                                      final_axes_number, 
                                                      GL_number, 
                                                      delais_TT_stop_del_axis, 
                                                      dTT_MS_cohort, 
                                                      TT_col_break,
                                                      user_parameters_completeness, 
                                                      user_dims_completeness)
    
    expected_axis_table = pandas.read_csv(full_full_expected_results_dir/'linear_second_axis_table.csv')
    expected_absolute_phen_table = pandas.read_csv(full_full_expected_results_dir/'linear_second_absolute_phen_table.csv')
    expected_relative_phen_table = pandas.read_csv(full_full_expected_results_dir/'linear_second_relative_phen_table.csv')
    expected_absolute_dim_table = pandas.read_csv(full_full_expected_results_dir/'linear_second_dim_table.csv')
    expected_parameters_table = pandas.read_csv(full_full_expected_results_dir/'linear_second_parameters_table.csv')
    expected_first_leaf_phen_table = pandas.read_csv(full_full_expected_results_dir/'linear_first_leaf_phen_table.csv')
    expected_HS_GL_SSI_dynamic_table = pandas.read_csv(full_full_expected_results_dir/'linear_HS_GL_SSI_dynamic_table.csv')
    expected_relative_dim_table = pandas.read_csv(full_full_expected_results_dir/'linear_relative_dim_table.csv')
    expected_tillering_dynamic_table = pandas.read_csv(full_full_expected_results_dir/'tillering_dynamic_table.csv')

    to_compare = {'linear_second_axis_table': (expected_axis_table, results[0]),
                 'linear_second_absolute_phen_table': (expected_absolute_phen_table, results[1]),
                 'linear_second_relative_phen_table': (expected_relative_phen_table, results[2]),
                 'linear_second_dim_table': (expected_absolute_dim_table, results[3]),
                 'linear_second_parameters_table': (expected_parameters_table, results[4]),
                 'linear_first_leaf_phen_table': (expected_first_leaf_phen_table, results[5]),
                 'linear_HS_GL_SSI_dynamic_table': (expected_HS_GL_SSI_dynamic_table, results[6]),
                 'linear_relative_dim_table': (expected_relative_dim_table, results[7]),
                 'tillering_dynamic_table': (expected_tillering_dynamic_table, results[8])}
    
    _check_results(to_compare, user_parameters_completeness, user_dims_completeness)


def _check_results(to_compare, user_parameters_completeness, user_dims_completeness):
    result_table_dir = fitting_results_directory.joinpath('%d_%d' % (user_parameters_completeness, user_dims_completeness))
    if not result_table_dir.exists():
        result_table_dir.mkdir()
    for key, value in to_compare.iteritems():
        expected_table = value[0]
        result_table = value[1]
        result_table_filepath = result_table_dir.joinpath(key + '.csv')
        result_table.to_csv(result_table_filepath, na_rep='NA', index=False)  
        print 'The results have been saved to %s' % result_table_filepath
        np.testing.assert_allclose(result_table.values, expected_table.values, relative_tolerance, absolute_tolerance)


# test the visualea node
pm = PackageManager()
pm.init(verbose=False)

def test_fit_adel_input_data():
    res = run(('alinea.adel.Tutorials','fit_adel_input_data'), {}, pm=pm)
    assert res == []
    

