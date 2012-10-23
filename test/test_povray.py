import pandas
import numpy as np
from openalea.core.path import path
from adel.povray import post_processing
    
def test_count_pixels():
    scene_path = path("data/test_povray/scene.bmp")
    scene_box_path = path("data/test_povray/scene_box.bmp")
    
    rgb_colors = [[0, 170, 0], 
                  [170, 170, 0], 
                  [255, 0, 0], 
                  [255, 170, 0], 
                  [255, 201, 5], 
                  [0, 0, 255], 
                  [255, 170, 0], 
                  [0, 0, 255], 
                  [255, 201, 5], 
                  [170, 255, 255], 
                  [170, 255, 255]]
    
    import tempfile
    povray_results_directory = path(tempfile.mkdtemp(suffix='_povray_results'))
    result_path = povray_results_directory.joinpath('scene_res.csv')
    normalized_scene_path = povray_results_directory.joinpath('normalized_scene.bmp')
    normalize = True
    post_processing.count_pixels(scene_path=scene_path, scene_box_path=scene_box_path, rgb_colors=rgb_colors, result_path=result_path, normalize=normalize, normalized_scene_path=normalized_scene_path)
    expected_result_array = np.array([['scene.bmp', 0.0, 28526.0, 375.0, 2253.0, 25.0, 83.0, 2959.0]], dtype=object)    
    result_array = pandas.read_csv(result_path).values
    np.testing.assert_equal(result_array, expected_result_array)
