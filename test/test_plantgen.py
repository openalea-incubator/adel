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

plants_number = 5
plants_density = 12
decide_child_axis_probabilities={'T0': 0.0, 'T1': 0.900, 'T2': 0.967, 'T3': 0.817, 'T4': 0.083}
decide_child_cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(decide_child_axis_probabilities)
MS_leaves_number_probabilities = {'10': 0.145, '11': 0.818, '12': 0.037, '13': 0.0, '14': 0.0}
ears_density = 25
GL_number = {1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}
delais_TT_stop_del_axis = 600
TT_col_break = 0.0
TT_col_N_phytomer_potential = {'MS': 1078.0, 'T1': 1148.0, 'T2': 1158.0, 'T3': 1168.0, 'T4': 1178.0}
number_of_ears = plants_number * ears_density / float(plants_density)

expected_results_dir = path('data/test_plantgen')
default_expected_results_dir = expected_results_dir.joinpath('default')
min_min_expected_results_dir = expected_results_dir.joinpath('MIN_MIN')
short_short_expected_results_dir = expected_results_dir.joinpath('SHORT_SHORT')
full_full_expected_results_dir = expected_results_dir.joinpath('FULL_FULL')

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
    expected_axeT = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    axeT_ = axeT.create_axeT_tmp(plants_number, decide_child_cohort_probabilities, MS_leaves_number_probabilities)
    test_table_filepath = default_results.joinpath('axeT_tmp.csv')
    axeT_.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    assert (axeT_['id_axis'] == expected_axeT['id_axis']).all() 
    del axeT_['id_axis']
    del expected_axeT['id_axis']
    np.testing.assert_allclose(axeT_.values, expected_axeT.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_dynT_tmp():
    axeT_ = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_dynT = pandas.read_csv(default_expected_results_dir/'dynT_tmp.csv')
    dynT_ = dynT.create_dynT_tmp(axeT_)
    test_table_filepath = default_results.joinpath('dynT_tmp.csv')
    dynT_.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    assert (dynT_['id_axis'] == expected_dynT['id_axis']).all()
    del dynT_['id_axis']
    del expected_dynT['id_axis']
    np.testing.assert_allclose(dynT_.values, expected_dynT.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_dimT_tmp():
    axeT_ = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_dimT = pandas.read_csv(default_expected_results_dir/'dimT_tmp.csv')
    dimT_ = dimT.create_dimT_tmp(axeT_)
    test_table_filepath = default_results.joinpath('dimT_tmp.csv')
    dimT_.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    assert (dimT_['id_axis'] == expected_dimT['id_axis']).all()
    del dimT_['id_axis']
    del expected_dimT['id_axis']
    np.testing.assert_allclose(dimT_.values, expected_dimT.values, relative_tolerance, absolute_tolerance)

        
@with_setup(reinit_random_state)
def test_create_dynT():
    dimT_tmp = pandas.read_csv(default_expected_results_dir/'dimT_tmp_merged.csv')
    dynT_tmp = pandas.read_csv(default_expected_results_dir/'dynT_tmp_merged.csv')
    decimal_elongated_internode_number = dynT.calculate_decimal_elongated_internode_number(dimT_tmp, dynT_tmp)
    expected_dynT = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    dynT_ = dynT.create_dynT(dynT_tmp, GL_number, decimal_elongated_internode_number)
    test_table_filepath = default_results.joinpath('dynT.csv')
    dynT_.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    assert (dynT_['id_axis'] == expected_dynT['id_axis']).all()
    del dynT_['id_axis']
    del expected_dynT['id_axis']
    np.testing.assert_allclose(dynT_.values, expected_dynT.values, relative_tolerance, absolute_tolerance)
    

@with_setup(reinit_random_state)
def test_create_phenT_tmp():
    dynT_ = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    axeT_tmp = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_phenT_tmp = pandas.read_csv(default_expected_results_dir/'phenT_tmp.csv')
    phenT_tmp = phenT.create_phenT_tmp(axeT_tmp, dynT_)
    test_table_filepath = default_results.joinpath('phenT_tmp.csv')
    phenT_tmp.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(phenT_tmp.values, expected_phenT_tmp.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_phenT_first():
    phenT_tmp = pandas.read_csv(default_expected_results_dir/'phenT_tmp.csv')
    expected_phenT_first = pandas.read_csv(default_expected_results_dir/'phenT_first.csv')
    phenT_first = phenT.create_phenT_first(phenT_tmp)
    test_table_filepath = default_results.joinpath('phenT_first.csv')
    phenT_first.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(phenT_first.values, expected_phenT_first.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_phenT():
    phenT_abs = pandas.read_csv(default_expected_results_dir/'phenT_abs.csv')
    expected_phenT = pandas.read_csv(default_expected_results_dir/'phenT.csv')
    phenT_first = pandas.read_csv(default_expected_results_dir/'phenT_first.csv')
    phenT_ = phenT.create_phenT(phenT_abs, phenT_first)
    test_table_filepath = default_results.joinpath('phenT.csv')
    phenT_.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(phenT_.values, expected_phenT.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_HS_GL_SSI_T():
    phenT_abs = pandas.read_csv(default_expected_results_dir/'phenT_abs.csv')
    axeT_tmp = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    dynT_ = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    expected_HS_GL_SSI_T = pandas.read_csv(default_expected_results_dir/'HS_GL_SSI_T.csv')
    HS_GL_SSI_T = phenT.create_HS_GL_SSI_T(phenT_abs, axeT_tmp, dynT_)
    test_table_filepath = default_results.joinpath('HS_GL_SSI_T.csv')
    HS_GL_SSI_T.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(HS_GL_SSI_T.values, expected_HS_GL_SSI_T.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_axeT():
    axeT_tmp = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_axeT = pandas.read_csv(default_expected_results_dir/'axeT.csv')
    phenT_first = pandas.read_csv(default_expected_results_dir/'phenT_first.csv')
    dynT_ = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    t1_most_frequent_MS = dynT_['t1'][dynT_.first_valid_index()]
    TT_regression_start = t1_most_frequent_MS + params.DELAIS_REG_MONT
    axeT_ = axeT.create_axeT(axeT_tmp, phenT_first, dynT_, TT_regression_start, TT_col_N_phytomer_potential['MS'], delais_TT_stop_del_axis, number_of_ears)
    test_table_filepath = default_results.joinpath('axeT.csv')
    axeT_.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    assert (axeT_['id_axis'] == expected_axeT['id_axis']).all()
    del axeT_['id_axis']
    del expected_axeT['id_axis']
    np.testing.assert_allclose(axeT_.values, expected_axeT.values, relative_tolerance, absolute_tolerance)
 

@with_setup(reinit_random_state)
def test_create_dimT_abs():
    axeT_ = pandas.read_csv(default_expected_results_dir/'axeT.csv')
    dimT_tmp = pandas.read_csv(default_expected_results_dir/'dimT_tmp_merged.csv')
    phenT_tmp = pandas.read_csv(default_expected_results_dir/'phenT_tmp.csv')
    dynT_ = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    expected_dimT = pandas.read_csv(default_expected_results_dir/'dimT_abs.csv')
    dimT_ = dimT.create_dimT_abs(axeT_, dimT_tmp, phenT_tmp, dynT_)
    test_table_filepath = default_results.joinpath('dimT_abs.csv')
    dimT_.to_csv(test_table_filepath, na_rep='NA', index=False)
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dimT_.values, expected_dimT.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_phenT_abs():
    phenT_tmp = pandas.read_csv(default_expected_results_dir/'phenT_tmp.csv')
    axeT_ = pandas.read_csv(default_expected_results_dir/'axeT.csv')
    dimT_abs = pandas.read_csv(default_expected_results_dir/'dimT_abs.csv')
    expected_phenT_abs = pandas.read_csv(default_expected_results_dir/'phenT_abs.csv')
    phenT_abs = phenT.create_phenT_abs(phenT_tmp, axeT_, dimT_abs)
    test_table_filepath = default_results.joinpath('phenT_abs.csv')
    phenT_abs.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(phenT_abs.values, expected_phenT_abs.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_tilleringT():
    axeT_ = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_tilleringT = pandas.read_csv(default_expected_results_dir/'tilleringT.csv')
    dynT_ = pandas.read_csv(default_expected_results_dir/'dynT.csv')
    phenT_first = pandas.read_csv(default_expected_results_dir/'phenT_first.csv')
    dynT_most_frequent_MS = dynT_.ix[0]
    id_cohort_most_frequent_MS = str(dynT_most_frequent_MS['id_cohort'])
    N_phytomer_potential_most_frequent_MS = str(dynT_most_frequent_MS['N_phytomer_potential'])
    id_phen_most_frequent_MS = int(''.join([id_cohort_most_frequent_MS, N_phytomer_potential_most_frequent_MS]))
    TT_start = phenT_first['TT_app_phytomer'][phenT_first[phenT_first['id_phen'] == id_phen_most_frequent_MS].index[0]]
    t1_most_frequent_MS = dynT_most_frequent_MS['t1']
    TT_regression_start = t1_most_frequent_MS + params.DELAIS_REG_MONT
    tilleringT = axeT.create_tilleringT(TT_start, TT_regression_start, TT_col_N_phytomer_potential['MS'], plants_number, plants_density, axeT_.index.size, ears_density)
    test_table_filepath = default_results.joinpath('tilleringT.csv')
    tilleringT.to_csv(test_table_filepath, na_rep='NA', index=False)
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(tilleringT.values, expected_tilleringT.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_create_cardinalityT():
    axeT_ = pandas.read_csv(default_expected_results_dir/'axeT_tmp.csv')
    expected_cardinalityT = pandas.read_csv(default_expected_results_dir/'cardinalityT.csv')
    (theoretical_cohort_cardinalities, 
     theoretical_axis_cardinalities) = tools.calculate_theoretical_cardinalities(
                                            plants_number, 
                                            decide_child_cohort_probabilities,
                                            decide_child_axis_probabilities,
                                            params.FIRST_CHILD_DELAY)
    cardinalityT = axeT.create_cardinalityT(theoretical_cohort_cardinalities, theoretical_axis_cardinalities, axeT_[['id_cohort', 'id_axis']])
    test_table_filepath = default_results.joinpath('cardinalityT.csv')
    cardinalityT.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    assert (cardinalityT['id_axis'] == expected_cardinalityT['id_axis']).all()
    del cardinalityT['id_axis']
    del expected_cardinalityT['id_axis']
    np.testing.assert_allclose(cardinalityT.values.astype(float), expected_cardinalityT.values, relative_tolerance, absolute_tolerance)
    

@with_setup(reinit_random_state)
def test_create_dimT():
    dimT_abs = pandas.read_csv(default_expected_results_dir/'dimT_abs.csv')
    expected_dimT = pandas.read_csv(default_expected_results_dir/'dimT.csv')
    dimT_ = dimT.create_dimT(dimT_abs)
    test_table_filepath = default_results.joinpath('dimT.csv')
    dimT_.to_csv(test_table_filepath, na_rep='NA', index=False)  
    print 'The results have been saved to %s' % test_table_filepath
    np.testing.assert_allclose(dimT_.values, expected_dimT.values, relative_tolerance, absolute_tolerance)


@with_setup(reinit_random_state)
def test_gen_adel_input_data_from_min():
    dynT_user = pandas.read_csv(min_min_expected_results_dir/'dynT_user.csv')
    dimT_user = pandas.read_csv(min_min_expected_results_dir/'dimT_user.csv')
    
    results = plantgen.gen_adel_input_data(dynT_user,
                                            dimT_user, 
                                            plants_number, 
                                            plants_density,
                                            decide_child_axis_probabilities, 
                                            MS_leaves_number_probabilities, 
                                            ears_density, 
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
    expected_cardinalityT = pandas.read_csv(min_min_expected_results_dir/'cardinalityT.csv')

    to_compare = {'axeT': (expected_axeT, results[0]),
                 'dimT': (expected_dimT, results[1]),
                 'phenT': (expected_phenT, results[2]),
                 'phenT_abs': (expected_phenT_abs, results[3]),
                 'dimT_abs': (expected_dimT_abs, results[4]),
                 'dynT': (expected_dynT, results[5]),
                 'phenT_first': (expected_phenT_first, results[6]),
                 'HS_GL_SSI_T': (expected_HS_GL_SSI_T, results[7]),
                 'tilleringT': (expected_tilleringT, results[8]),
                 'cardinalityT': (expected_cardinalityT, results[9])}
    
    _check_results(to_compare, plantgen.DataCompleteness.MIN, plantgen.DataCompleteness.MIN)

         
@with_setup(reinit_random_state)
def test_gen_adel_input_data_from_short():
    dynT_user = pandas.read_csv(short_short_expected_results_dir/'dynT_user.csv')
    dimT_user = pandas.read_csv(short_short_expected_results_dir/'dimT_user.csv')
    results = plantgen.gen_adel_input_data(dynT_user,
                                          dimT_user, 
                                          plants_number, 
                                          plants_density,
                                          decide_child_axis_probabilities, 
                                          MS_leaves_number_probabilities, 
                                          ears_density, 
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
    expected_cardinalityT = pandas.read_csv(short_short_expected_results_dir/'cardinalityT.csv')

    to_compare = {'axeT': (expected_axeT, results[0]),
                  'dimT': (expected_dimT, results[1]),
                  'phenT': (expected_phenT, results[2]),
                  'phenT_abs': (expected_phenT_abs, results[3]),
                  'dimT_abs': (expected_dimT_abs, results[4]),
                  'dynT': (expected_dynT, results[5]),
                  'phenT_first': (expected_phenT_first, results[6]),
                  'HS_GL_SSI_T': (expected_HS_GL_SSI_T, results[7]),
                  'tilleringT': (expected_tilleringT, results[8]),
                  'cardinalityT': (expected_cardinalityT, results[9])}
    
    _check_results(to_compare, plantgen.DataCompleteness.SHORT, plantgen.DataCompleteness.SHORT)


@with_setup(reinit_random_state)
def test_gen_adel_input_data_from_full():
    dynT_user = pandas.read_csv(full_full_expected_results_dir/'dynT_user.csv')
    dimT_user = pandas.read_csv(full_full_expected_results_dir/'dimT_user.csv')
    results = plantgen.gen_adel_input_data(dynT_user,
                                         dimT_user,
                                         plants_number, 
                                         plants_density,
                                         decide_child_axis_probabilities, 
                                         MS_leaves_number_probabilities, 
                                         ears_density, 
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
    expected_cardinalityT = pandas.read_csv(full_full_expected_results_dir/'cardinalityT.csv')

    to_compare = {'axeT': (expected_axeT, results[0]),
                 'dimT': (expected_dimT, results[1]),
                 'phenT': (expected_phenT, results[2]),
                 'phenT_abs': (expected_phenT_abs, results[3]),
                 'dimT_abs': (expected_dimT_abs, results[4]),
                 'dynT': (expected_dynT, results[5]),
                 'phenT_first': (expected_phenT_first, results[6]),
                 'HS_GL_SSI_T': (expected_HS_GL_SSI_T, results[7]),
                 'tilleringT': (expected_tilleringT, results[8]),
                 'cardinalityT': (expected_cardinalityT, results[9])}
    
    _check_results(to_compare, plantgen.DataCompleteness.FULL, plantgen.DataCompleteness.FULL)


def _check_results(to_compare, dynT_user_completeness, dimT_user_completeness):
    result_table_dir = tmp_results_directory.joinpath('%s_%s' % (dynT_user_completeness, 
                                                                 dimT_user_completeness))
    if not result_table_dir.exists():
        result_table_dir.mkdir()
    for key, value in to_compare.iteritems():
        expected_table = value[0]
        result_table = value[1]
        result_table_filepath = result_table_dir.joinpath(key + '.csv')
        result_table.to_csv(result_table_filepath, na_rep='NA', index=False)  
        print 'The results have been saved to %s' % result_table_filepath
        if 'id_axis' in expected_table:
            assert (result_table['id_axis'] == expected_table['id_axis']).all()
            del result_table['id_axis']
            del expected_table['id_axis']
        np.testing.assert_allclose(result_table.values, expected_table.values, relative_tolerance, absolute_tolerance)


def test_plantgen():
    # test the visualea node
    pm = PackageManager()
    pm.init(verbose=False)
    res = run(('alinea.adel.Tutorials', 'plantgen'), {}, pm=pm)
    assert res == []
    
