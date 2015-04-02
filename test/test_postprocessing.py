import tempfile

import numpy as np
import pandas as pd
from openalea.core.path import path

from alinea.adel import postprocessing as pp

DATA_DIRPATH = path('data/test_postprocessing')

INPUTS_DIRPATH = DATA_DIRPATH/'inputs'

OUTPUTS_DIRPATH = DATA_DIRPATH/'outputs'

RELATIVE_TOLERANCE = 10e-3
ABSOLUTE_TOLERANCE = 10e-3


def test_aggregate_adel_output():
    adel_output_df = pd.read_csv(INPUTS_DIRPATH/'adel_output.csv')
    adel_output_aggregated_df = pp.aggregate_adel_output(adel_output_df, by=['TT', 'plant', 'axe_id'])
    adel_output_aggregated_df.to_csv(OUTPUTS_DIRPATH/'actual_adel_aggregated_output.csv', index=False, na_rep='NA')
    desired_adel_output_aggregated_df = pd.read_csv(OUTPUTS_DIRPATH/'desired_adel_aggregated_output.csv')
    adel_output_aggregated_df = adel_output_aggregated_df.select_dtypes(include=[np.number])
    desired_adel_output_aggregated_df = desired_adel_output_aggregated_df.select_dtypes(include=[np.number])
    np.testing.assert_allclose(adel_output_aggregated_df.values, desired_adel_output_aggregated_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_phenology():
    adel_output_df = pd.read_csv(INPUTS_DIRPATH/'adel_output.csv')
    phenology_df = pp.phenology(adel_output_df)
    phenology_df.to_csv(OUTPUTS_DIRPATH/'actual_phenology.csv', index=False, na_rep='NA')
    desired_phenology_df = pd.read_csv(OUTPUTS_DIRPATH/'desired_phenology.csv')
    phenology_df = phenology_df.select_dtypes(include=[np.number])
    desired_phenology_df = desired_phenology_df.select_dtypes(include=[np.number])
    np.testing.assert_allclose(phenology_df.values, desired_phenology_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)
 
 
def test_axis_statistics():
    adel_output_df = pd.read_csv(INPUTS_DIRPATH/'adel_output.csv')
    axis_statistics_df = pp.axis_statistics(adel_output_df, domain_area=1)
    axis_statistics_df.to_csv(OUTPUTS_DIRPATH/'actual_axis_statistics.csv', index=False, na_rep='NA')
    desired_axis_statistics_df = pd.read_csv(OUTPUTS_DIRPATH/'desired_axis_statistics.csv')
    axis_statistics_df = axis_statistics_df.select_dtypes(include=[np.number])
    desired_axis_statistics_df = desired_axis_statistics_df.select_dtypes(include=[np.number])
    np.testing.assert_allclose(axis_statistics_df.values, desired_axis_statistics_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)
 
 
def test_plot_statistics():
    axis_statistics_df = pd.read_csv(INPUTS_DIRPATH/'axis_statistics.csv')
    plot_statistics_df = pp.plot_statistics(axis_statistics_df, plant_number=1, domain_area=1)
    plot_statistics_df.to_csv(OUTPUTS_DIRPATH/'actual_plot_statistics.csv', index=False, na_rep='NA')
    desired_plot_statistics_df = pd.read_csv(OUTPUTS_DIRPATH/'desired_plot_statistics.csv')
    plot_statistics_df = plot_statistics_df.select_dtypes(include=[np.number])
    desired_plot_statistics_df = desired_plot_statistics_df.select_dtypes(include=[np.number])
    np.testing.assert_allclose(plot_statistics_df.values, desired_plot_statistics_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)

