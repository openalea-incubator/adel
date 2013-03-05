import tempfile

import numpy as np
import pandas
from openalea.core.path import path

from alinea.adel.stand.stand import post_processing

relative_tolerance = 10e-3
absolute_tolerance = 10e-3
    
def test_post_processing():
    adel_output_path = path('data/test_stand/adel_output.csv')
    
    import tempfile
    global_postprocessing_file_path = path(tempfile.mktemp(suffix='.csv'))
    peraxis_postprocessing_file_path = path(tempfile.mktemp(suffix='.csv'))
    
    # launch the postprocessing
    (global_postprocessing_file_path, 
     peraxis_postprocessing_file_path,
     intermediate_results_file_path)= \
        post_processing(adel_output_path, 
                        1, 
                        0.0046082949, 
                        global_postprocessing_file_path,
                        peraxis_postprocessing_file_path)
        
    expected_results_dir_path = path('data/test_stand')
    # check intermediate results
    intermediate_results_array = \
        pandas.read_csv(intermediate_results_file_path).values
    expected_intermediate_results_file_path = \
        expected_results_dir_path/'intermediate.csv'
    expected_intermediate_results_array = \
        pandas.read_csv(expected_intermediate_results_file_path).values
    np.testing.assert_allclose(intermediate_results_array, 
                               expected_intermediate_results_array,
                               relative_tolerance, 
                               absolute_tolerance)
    
    for (expected_results, results_to_test) \
        in [(expected_results_dir_path/'global_postprocessing.csv', global_postprocessing_file_path),
            (expected_results_dir_path/'peraxis_postprocessing.csv', peraxis_postprocessing_file_path)]:
        # check postprocessing results
        postprocessing_results_df = pandas.read_csv(results_to_test)
        postprocessing_results_df = \
            postprocessing_results_df.take(
                range(postprocessing_results_df.columns.size)[1:], 
                1)
        postprocessing_results_array = postprocessing_results_df.values
        expected_postprocessing_results_df = pandas.read_csv(expected_results)
        expected_postprocessing_results_array = expected_postprocessing_results_df.values
        np.testing.assert_allclose(postprocessing_results_array, 
                                   expected_postprocessing_results_array,
                                   relative_tolerance, 
                                   absolute_tolerance)
    
if __name__ == '__main__':
    test_post_processing()
    
