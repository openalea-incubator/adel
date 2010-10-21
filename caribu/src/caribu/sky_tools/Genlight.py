import os

class Genlight(object):
    """  Doc... """ 

    def __init__(self):
        pass


    def __call__(self, SkyObj,filename):
        """ SkyObj = table contenant pour chaque ligne [I, Vx,Vy,Vz]  """
        if filename=='sky.light':#nom de fichier par defaut, sauve dans un dossier temporaire
            d=os.tempnam()
            os.mkdir(d)
            pathname = os.path.join(d,filename)
        else:# dans le dossier courant avec nom de fichier personnalise - a mettre a jour
            pathname = filename

        f = file(pathname, 'w')
        for i in range(len(SkyObj.sky)):
            f.write(str(SkyObj.sky[i][0])+' '+str(SkyObj.sky[i][1])+' '+str(SkyObj.sky[i][2])+' '+str(SkyObj.sky[i][3])+'\n')

        f.close()

        return pathname



