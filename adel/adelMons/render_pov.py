import os
from os.path import join
from os.path import split
import win32api
import time

path_ini = os.getcwd()

class render_pov(object):
    """  render a pov file and generates a bmp file . !! Pov-ray has to be launch !! """ 

    def __init__(self):
        pass


    def __call__(self, pov_file, povdir, wait):

        dir = pov_file[:-len(pov_file.split('\\')[-1])-1]#path du fichier .pov d'origine
        name = pov_file.split('\\')[-1].split('.')[0]

         
        #os.system('copy '+win32api.GetShortPathName(pov_file)+' '+win32api.GetShortPathName(povdir))  #copie du fichier pov dans le bin de povray pas necessaire : avec short paths peut executer diretcement
        os.chdir(povdir)
        os.system('pvengine.exe /RENDER '+join(win32api.GetShortPathName(dir),pov_file.split('\\')[-1]))#short name pour le path, mais pas le nom de fichier pour conserver noms'd'image complets
        time.sleep(wait) #temps d'attente en secondes entre 2 simulations
        os.chdir(path_ini)
        
        return join(dir, name+'.bmp')

