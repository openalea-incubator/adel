import Sky

class GenSky(object):
    """  Doc... """ 

    def __init__(self):
        pass


    def __call__(self, Rd, Tsky, Nazi, Nzen):
        """ renvoie liste d'intensites associees aux secteurs zenitaux et azimutaux """
        SkyRd=Sky.Sky(Nazi, Nzen)
        SkyRd.set_Rd(Rd,Tsky)

        return SkyRd#(SkyRd.sky,)#pour renvoyer liste en entier et pas seulement le 1er element: lourdeur
        
