from openalea.plantgl.all import Viewer
import os.path

def save_image(scene, image_name='%s/img%d.%s', directory='.', index=0, ext='png'):
    '''
    Save an image of a scene in a specific directory

    Parameters
    ----------

        - scene: a PlantGL scene
        - image_name: a string template 
            The format of the string is dir/img5.png
        - directory (optional: ".") the directory where the images are written
        - index: the index of the image
        - ext : the image format
    '''

    filename = image_name%(directory, index, ext)
    Viewer.frameGL.saveImage(filename)
    return scene,
