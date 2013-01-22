import tempfile

import numpy as np
import pandas
from openalea.core.path import path

from alinea.adel.stand.stand import post_processing
    
def test_post_processing():
    adel_output_path = path("data/test_stand/adel_output.csv")
    
    import tempfile
    postprocessing_results_file_path = path(tempfile.mktemp(suffix='.csv'))
    expected_postprocessing_results_file_path = \
        path("data/test_stand/postprocessing.csv")
    intermediate_results_file_path = path(tempfile.mktemp(suffix='.csv'))
    expected_intermediate_results_file_path = \
        path("data/test_stand/intermediate.csv")
        
    (postprocessing_results_file_path, 
     intermediate_results_file_path)= \
        post_processing(adel_output_path, 
                        218, 
                        1.00909090909, 
                        postprocessing_results_file_path)
    intermediate_results_array = \
        pandas.read_csv(intermediate_results_file_path).values
    expected_intermediate_results_array = \
        pandas.read_csv(expected_intermediate_results_file_path).values
    np.testing.assert_allclose(intermediate_results_array, 
                            expected_intermediate_results_array)
    postprocessing_results_df = pandas.read_csv(postprocessing_results_file_path)
    postprocessing_results_df = \
        postprocessing_results_df.take(
            range(postprocessing_results_df.columns.size)[1:], 
            1)
    postprocessing_results_array = postprocessing_results_df.values
    expected_postprocessing_results_df = pandas.read_csv(expected_postprocessing_results_file_path)
    expected_postprocessing_results_array = expected_postprocessing_results_df.values
    np.testing.assert_allclose(postprocessing_results_array, expected_postprocessing_results_array)

