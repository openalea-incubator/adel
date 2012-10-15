import pandas
import numpy as np
from PIL import Image
from openalea.core.path import path

WHITE = (255, 255, 255)

def count_pixels(scene_path='', scene_box_path='', colors=[], result_path=''):
    '''Counts the pixels of each color.
    :Parameters:
        - `scene_path` : the path of the bmp file which represents the scene 
        image.
        - `scene_box_path` : the path of the bmp file which represents the area 
        where the pixels have to be counted.
        - `colors` : the colors to consider when counting the pixels. Each 
        color is described by a 3-tuple: (R, G, B). 
        - `result_path` : the path of csv file where the result is added. 
        The first column contains the base name of the bmp file which 
        represents the scene image (i.e. scene_path). Other 
        columns are the list of colors to consider.  
    :Types:
        - `scene_path` : str     
        - `scene_box_path` : str     
        - `colors` : list of 3-tuple of int  
        - `result_path` : str     
        
    :return: the path of the csv file where the result is added. See the 
    documentation of the input arguments for information about the structure.
    :rtype: str
    ''' 
    assert scene_path != '' and scene_box_path != '' and len(colors) != 0 and result_path != ''
    scene_path = path(scene_path)
    scene_box_path = path(scene_box_path)
    result_path = path(result_path)
    
    scene_image = Image.open(scene_path)
    scene_box_image = Image.open(scene_box_path)
    
    scene_array = scene_image.load()
    scene_box_array = scene_box_image.load()
    width, height = scene_image.size
    
    colors_tuple = set([tuple(color) for color in colors])
    pixels_init_number = np.zeros_like(range(len(colors_tuple)))
    pixels_numbers = dict(zip(colors_tuple, pixels_init_number))
    
    for x in range(width):
        for y in range(height):
            if scene_box_array[x, y] == WHITE:
                color = scene_array[x, y]
                if color in colors_tuple:
                    pixels_numbers[color] += 1
    
    if result_path.exists():
        result_df = pandas.read_csv(result_path)
    else:
        columns = ['Filename'] + [str(color) for color in colors_tuple]
        result_df = pandas.DataFrame(columns=columns)
    
    data = [[scene_path.basename()] + [pixels_numbers[color] for color in colors_tuple]]
    new_index = [result_df.index.size]
    new_df = pandas.DataFrame(data, index=new_index, columns=result_df.columns)
    
    result_df = pandas.concat([result_df, new_df])
    result_df.to_csv(result_path, na_rep='NA', index=False)
    
    return str(result_path)



