import random

from adel.plantgen import axeT, dimT, phenT, \
    plantgen, dynT, params, tools
import numpy as np
import pandas
from openalea.core.path import path
from openalea.core.alea import *
from nose.tools import with_setup

random.seed(1234)

initial_random_state = random.getstate()

plant_number = 5
decide_child_cohort_probabilities = {'3': 0.0, '4': 0.900, '5': 0.967, '6': 0.817, '7': 0.083}
MS_leaves_number_probabilities = {'10': 0.145, '11': 0.818, '12': 0.037, '13': 0.0, '14': 0.0}
TT_bolting = 500.0
final_axes_density = 15
GL_number = {1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}
delais_TT_stop_del_axis = 600
TT_col_nff = {'1': 1078.0, '4': 1148.0, '5': 1158.0, '6': 1168.0, '7': 1178.0}

expected_results_dir = path('data/test_plantgen')
default_expected_results_dir = expected_results_dir.joinpath('default')
min_min_expected_results_dir = expected_results_dir.joinpath('min_min')
short_short_expected_results_dir = expected_results_dir.joinpath('short_short')
full_full_expected_results_dir = expected_results_dir.joinpath('full_full')

completeness_mapping = {1: 'min', 2: 'short', 3: 'full'}

import tempfile
tmp_results_directory = path(tempfile.mkdtemp(suffix='_plantgen_results'))
default_results = tmp_results_directory.joinpath('default')
if not default_results.exists():
    default_results.mkdir()
relative_tolerance = 10e-3
absolute_tolerance = 10e-3


def reinit_random_state():
    random.setstate(initial_random_state)


@with_setup(reinit_random_state)
def test_create_axeT_tmp():
    expected_axeT_dataframe = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    axeT_dataframe = axeT.create_axeT_tmp(plant_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities)
    test_table_filepath = default_results.joinpath('axeT_tmp.csv')
    axeT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    assert (axeT_dataframe['id_axis'] == expected_axeT_dataframe['id_axis']).all() 
    del axeT_dataframe['id_axis']
    del expected_axeT_dataframe['id_axis']
    np.testing.assert_allclose(axeT_dataframe.values, expected_axeT_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_dynT_tmp():
    axeT_dataframe = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_dynT_dataframe = pandas.read_csv(default_expected_results_dir/'dynT_tmp.csv')
    dynT_dataframe = dynT.create_dynT_tmp(axeT_dataframe['id_phen'].tolist())
    test_table_filepath = default_results.joinpath('dynT_tmp.csv')
    dynT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dynT_dataframe.values, expected_dynT_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_dimT_tmp():
    dynT_dataframe = pandas.read_csv(default_expected_results_dir/'dynT_tmp.csv')
    expected_dimT_dataframe = pandas.read_csv(default_expected_results_dir/'dimT_tmp.csv')
    dimT_dataframe = dimT.create_dimT_tmp(dynT_dataframe)
    test_table_filepath = default_results.joinpath('dimT_tmp.csv')
    dimT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dimT_dataframe.values, expected_dimT_dataframe.values, relative_tolerance, absolute_tolerance)

        
@with_setup(reinit_random_state)
def test_create_dynT():
    dynT_user_dataframe = pandas.read_csv(default_expected_results_dir/'dynT_user.csv')
    dimT_user_dataframe = pandas.read_csv(default_expected_results_dir/'dimT_user.csv')
    decimal_elongated_internode_number = dynT.calculate_decimal_elongated_internode_number(dimT_user_dataframe)
    expected_dynT_dataframe = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    dynT_dataframe = dynT.create_dynT(dynT_user_dataframe, dimT_user_dataframe, GL_number, decimal_elongated_internode_number)
    test_table_filepath = default_results.joinpath('dynT.csv')
    dynT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dynT_dataframe.values, expected_dynT_dataframe.values, relative_tolerance, absolute_tolerance)
    

@with_setup(reinit_random_state)
def test_create_phenT_abs():
    dynT_dataframe = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    dimT_user_dataframe = pandas.read_csv(default_expected_results_dir/'dimT_user.csv')
    decimal_elongated_internode_number = dynT.calculate_decimal_elongated_internode_number(dimT_user_dataframe)
    expected_phenT_abs_dataframe = pandas.read_csv(default_expected_results_dir/'phenT_abs.csv')
    phenT_abs_dataframe = phenT.create_phenT_abs(dynT_dataframe, decimal_elongated_internode_number)
    test_table_filepath = default_results.joinpath('phenT_abs.csv')
    phenT_abs_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(phenT_abs_dataframe.values, expected_phenT_abs_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_phenT_first():
    phenT_abs_dataframe = pandas.read_csv(default_expected_results_dir/'phenT_abs.csv')
    expected_phenT_first_dataframe = pandas.read_csv(default_expected_results_dir/'phenT_first.csv')
    phenT_first_dataframe = phenT.create_phenT_first(phenT_abs_dataframe)
    test_table_filepath = default_results.joinpath('phenT_first.csv')
    phenT_first_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(phenT_first_dataframe.values, expected_phenT_first_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_phenT():
    phenT_abs_dataframe = pandas.read_csv(default_expected_results_dir/'phenT_abs.csv')
    expected_phenT_dataframe = pandas.read_csv(default_expected_results_dir/'phenT.csv')
    phenT_first_dataframe = pandas.read_csv(default_expected_results_dir/'phenT_first.csv')
    phenT_dataframe = phenT.create_phenT(phenT_abs_dataframe, phenT_first_dataframe)
    test_table_filepath = default_results.joinpath('phenT.csv')
    phenT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(phenT_dataframe.values, expected_phenT_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_HS_GL_SSI_T():
    dynT_dataframe = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    expected_HS_GL_SSI_T_dataframe = pandas.read_csv(default_expected_results_dir/'HS_GL_SSI_T.csv')
    HS_GL_SSI_T_dataframe = phenT.create_HS_GL_SSI_T(dynT_dataframe)
    test_table_filepath = default_results.joinpath('HS_GL_SSI_T.csv')
    HS_GL_SSI_T_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(HS_GL_SSI_T_dataframe.values, expected_HS_GL_SSI_T_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_axeT():
    axeT_tmp_dataframe = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_axeT_dataframe = pandas.read_csv(default_expected_results_dir/'axeT.csv')
    phenT_first_dataframe = pandas.read_csv(default_expected_results_dir/'phenT_first.csv')
    dynT_dataframe = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    axeT_dataframe = axeT.create_axeT(axeT_tmp_dataframe, phenT_first_dataframe, dynT_dataframe, TT_bolting, TT_col_nff['1'], delais_TT_stop_del_axis, final_axes_density)
    test_table_filepath = default_results.joinpath('axeT.csv')
    axeT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    assert (axeT_dataframe['id_axis'] == expected_axeT_dataframe['id_axis']).all()
    axeT_dataframe = axeT_dataframe.select(lambda x: x not in ['id_axis'], 1)
    expected_axeT_dataframe = expected_axeT_dataframe.select(lambda x: x not in ['id_axis'], 1)
    np.testing.assert_allclose(axeT_dataframe.values, expected_axeT_dataframe.values, relative_tolerance, absolute_tolerance)
 

@with_setup(reinit_random_state)
def test_create_dimT_abs():
    axeT_dataframe = pandas.read_csv(default_expected_results_dir/'axeT.csv')
    dimT_user_dataframe = pandas.read_csv(default_expected_results_dir/'dimT_user.csv')
    phenT_abs_dataframe = pandas.read_csv(default_expected_results_dir/'phenT_abs.csv')
    dynT_dataframe = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    expected_dimT_dataframe = pandas.read_csv(default_expected_results_dir/'dimT_abs.csv')
    dimT_dataframe = dimT.create_dimT_abs(axeT_dataframe, dimT_user_dataframe, phenT_abs_dataframe, dynT_dataframe)
    test_table_filepath = default_results.joinpath('dimT_abs.csv')
    dimT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dimT_dataframe.values, expected_dimT_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_tilleringT():
    axeT_dataframe = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_tilleringT_dataframe = pandas.read_csv(default_expected_results_dir/'tilleringT.csv')
    tilleringT_dataframe = axeT.create_tilleringT(0, TT_bolting, TT_col_nff['1'], plant_number, axeT_dataframe, final_axes_density)
    test_table_filepath = default_results.joinpath('tilleringT.csv')
    tilleringT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(tilleringT_dataframe.values, expected_tilleringT_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_cohortT():
    axeT_dataframe = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_cohortT_dataframe = pandas.read_csv(default_expected_results_dir/'cohortT.csv')
    theoretical_cohorts_cardinalities = tools.calculate_theoretical_cohorts_cardinalities(plant_number, 
                                                                                          decide_child_cohort_probabilities,
                                                                                          params.FIRST_CHILD_DELAY)
    cohortT_dataframe = axeT.create_cohortT(theoretical_cohorts_cardinalities, axeT_dataframe['id_cohort_axis'])
    test_table_filepath = default_results.joinpath('cohortT.csv')
    cohortT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(cohortT_dataframe.values, expected_cohortT_dataframe.values, relative_tolerance, absolute_tolerance)
    

@with_setup(reinit_random_state)
def test_create_dimT():
    dimT_abs_dataframe = pandas.read_csv(default_expected_results_dir/'dimT_abs.csv')
    expected_dimT_dataframe = pandas.read_csv(default_expected_results_dir/'dimT.csv')
    dimT_dataframe = dimT.create_dimT(dimT_abs_dataframe)
    test_table_filepath = default_results.joinpath('dimT.csv')
    dimT_dataframe.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dimT_dataframe.values, expected_dimT_dataframe.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_gen_adel_input_data_from_min():
    dynT_user_completeness = plantgen.DataCompleteness.MIN
    dimT_user_completeness = plantgen.DataCompleteness.MIN
    TT_col_break = 0.0
    dynT_user = {'a_cohort': 0.0102, 
                 'TT_col_0': -0.771289027, 
                 'n0': 4.871559739,
                 'n1': 3.24283148,
                 'n2': 5.8}
    dimT_user = pandas.read_csv(min_min_expected_results_dir/'dimT_user.csv')
    results = plantgen.gen_adel_input_data_from_min(dynT_user,
                                                TT_col_nff,
                                                dimT_user, 
                                                plant_number, 
                                                decide_child_cohort_probabilities, 
                                                MS_leaves_number_probabilities, 
                                                TT_bolting, 
                                                final_axes_density, 
                                                GL_number, 
                                                delais_TT_stop_del_axis, 
                                                TT_col_break)
    
    expected_axeT = pandas.read_csv(min_min_expected_results_dir/'axeT.csv')
    expected_phenT_abs = pandas.read_csv(min_min_expected_results_dir/'phenT_abs.csv')
    expected_phenT = pandas.read_csv(min_min_expected_results_dir/'phenT.csv')
    expected_dimT_abs = pandas.read_csv(min_min_expected_results_dir/'dimT_abs.csv')
    expected_dynT = pandas.read_csv(min_min_expected_results_dir/'dynT.csv')
    expected_phenT_first = pandas.read_csv(min_min_expected_results_dir/'phenT_first.csv')
    expected_HS_GL_SSI_T = pandas.read_csv(min_min_expected_results_dir/'HS_GL_SSI_T.csv')
    expected_dimT = pandas.read_csv(min_min_expected_results_dir/'dimT.csv')
    expected_tilleringT = pandas.read_csv(min_min_expected_results_dir/'tilleringT.csv')
    expected_cohortT = pandas.read_csv(min_min_expected_results_dir/'cohortT.csv')

    to_compare = {'axeT': (expected_axeT, results[0]),
                 'dimT': (expected_dimT, results[1]),
                 'phenT': (expected_phenT, results[2]),
                 'phenT_abs': (expected_phenT_abs, results[3]),
                 'dimT_abs': (expected_dimT_abs, results[4]),
                 'dynT': (expected_dynT, results[5]),
                 'phenT_first': (expected_phenT_first, results[6]),
                 'HS_GL_SSI_T': (expected_HS_GL_SSI_T, results[7]),
                 'tilleringT': (expected_tilleringT, results[8]),
                 'cohortT': (expected_cohortT, results[9])}
    
    _check_results(to_compare, dynT_user_completeness, dimT_user_completeness)

         
@with_setup(reinit_random_state)
def test_gen_adel_input_data_from_short():
    dynT_user_completeness = plantgen.DataCompleteness.SHORT
    dimT_user_completeness = plantgen.DataCompleteness.SHORT
    TT_col_break = 0.0
    dynT_user = pandas.read_csv(short_short_expected_results_dir/'dynT_user.csv')
    dimT_user = pandas.read_csv(short_short_expected_results_dir/'dimT_user.csv')
    results = plantgen.gen_adel_input_data_from_short(dynT_user,
                                                  dimT_user, 
                                                  plant_number, 
                                                  decide_child_cohort_probabilities, 
                                                  MS_leaves_number_probabilities, 
                                                  TT_bolting, 
                                                  final_axes_density, 
                                                  GL_number, 
                                                  delais_TT_stop_del_axis, 
                                                  TT_col_break)
    
    expected_axeT = pandas.read_csv(short_short_expected_results_dir/'axeT.csv')
    expected_phenT_abs = pandas.read_csv(short_short_expected_results_dir/'phenT_abs.csv')
    expected_phenT = pandas.read_csv(short_short_expected_results_dir/'phenT.csv')
    expected_dimT_abs = pandas.read_csv(short_short_expected_results_dir/'dimT_abs.csv')
    expected_dynT = pandas.read_csv(short_short_expected_results_dir/'dynT.csv')
    expected_phenT_first = pandas.read_csv(short_short_expected_results_dir/'phenT_first.csv')
    expected_HS_GL_SSI_T = pandas.read_csv(short_short_expected_results_dir/'HS_GL_SSI_T.csv')
    expected_dimT = pandas.read_csv(short_short_expected_results_dir/'dimT.csv')
    expected_tilleringT = pandas.read_csv(short_short_expected_results_dir/'tilleringT.csv')
    expected_cohortT = pandas.read_csv(short_short_expected_results_dir/'cohortT.csv')

    to_compare = {'axeT': (expected_axeT, results[0]),
                 'dimT': (expected_dimT, results[1]),
                 'phenT': (expected_phenT, results[2]),
                 'phenT_abs': (expected_phenT_abs, results[3]),
                 'dimT_abs': (expected_dimT_abs, results[4]),
                 'dynT': (expected_dynT, results[5]),
                 'phenT_first': (expected_phenT_first, results[6]),
                 'HS_GL_SSI_T': (expected_HS_GL_SSI_T, results[7]),
                 'tilleringT': (expected_tilleringT, results[8]),
                 'cohortT': (expected_cohortT, results[9])}
    
    _check_results(to_compare, dynT_user_completeness, dimT_user_completeness)


@with_setup(reinit_random_state)
def test_gen_adel_input_data_from_full():
    dynT_user_completeness = plantgen.DataCompleteness.FULL
    dimT_user_completeness = plantgen.DataCompleteness.FULL
    TT_col_break = 0.0
    dynT_user = pandas.read_csv(full_full_expected_results_dir/'dynT_user.csv')
    dimT_user = pandas.read_csv(full_full_expected_results_dir/'dimT_user.csv')
    results = plantgen.gen_adel_input_data_from_full(dynT_user,
                                                 dimT_user,
                                                 plant_number, 
                                                 decide_child_cohort_probabilities, 
                                                 MS_leaves_number_probabilities, 
                                                 TT_bolting, 
                                                 final_axes_density, 
                                                 GL_number, 
                                                 delais_TT_stop_del_axis, 
                                                 TT_col_break)
    
    expected_axeT = pandas.read_csv(full_full_expected_results_dir/'axeT.csv')
    expected_phenT_abs = pandas.read_csv(full_full_expected_results_dir/'phenT_abs.csv')
    expected_phenT = pandas.read_csv(full_full_expected_results_dir/'phenT.csv')
    expected_dimT_abs = pandas.read_csv(full_full_expected_results_dir/'dimT_abs.csv')
    expected_dynT = pandas.read_csv(full_full_expected_results_dir/'dynT.csv')
    expected_phenT_first = pandas.read_csv(full_full_expected_results_dir/'phenT_first.csv')
    expected_HS_GL_SSI_T = pandas.read_csv(full_full_expected_results_dir/'HS_GL_SSI_T.csv')
    expected_dimT = pandas.read_csv(full_full_expected_results_dir/'dimT.csv')
    expected_tilleringT = pandas.read_csv(full_full_expected_results_dir/'tilleringT.csv')
    expected_cohortT = pandas.read_csv(full_full_expected_results_dir/'cohortT.csv')

    to_compare = {'axeT': (expected_axeT, results[0]),
                 'dimT': (expected_dimT, results[1]),
                 'phenT': (expected_phenT, results[2]),
                 'phenT_abs': (expected_phenT_abs, results[3]),
                 'dimT_abs': (expected_dimT_abs, results[4]),
                 'dynT': (expected_dynT, results[5]),
                 'phenT_first': (expected_phenT_first, results[6]),
                 'HS_GL_SSI_T': (expected_HS_GL_SSI_T, results[7]),
                 'tilleringT': (expected_tilleringT, results[8]),
                 'cohortT': (expected_cohortT, results[9])}
    
    _check_results(to_compare, dynT_user_completeness, dimT_user_completeness)


def _check_results(to_compare, dynT_user_completeness, dimT_user_completeness):
    result_table_dir = tmp_results_directory.joinpath('%s_%s' % (completeness_mapping[dynT_user_completeness], 
                                                                 completeness_mapping[dimT_user_completeness]))
    if not result_table_dir.exists():
        result_table_dir.mkdir()
    for key, value in to_compare.iteritems():
        expected_table = value[0]
        result_table = value[1]
        result_table_filepath = result_table_dir.joinpath(key + '.csv')
        result_table.to_csv(result_table_filepath, na_rep='NA', index=False)  
        print 'The results have been saved to %s' % result_table_filepath
        if key == 'axeT':
            assert (result_table['id_axis'] == expected_table['id_axis']).all()
            result_table = result_table.select(lambda x: x not in ['id_axis'], 1)
            expected_table = expected_table.select(lambda x: x not in ['id_axis'], 1)
        np.testing.assert_allclose(result_table.values, expected_table.values, relative_tolerance, absolute_tolerance)


def test_plantgen_MIN():
    # test the visualea node
    pm = PackageManager()
    pm.init(verbose=False)
    res = run(('alinea.adel.Tutorials', 'plantgen_MIN'), {}, pm=pm)
    assert res == []
    
