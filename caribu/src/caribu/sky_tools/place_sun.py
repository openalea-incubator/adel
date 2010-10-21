import Sky

class place_sun(object):
    """  Doc... """ 

    def __init__(self):
        pass


    def __call__(self, nbp, nbt, sun):
        s = Sky.Sky(nbp, nbt)
        s.set_Rsun(sun)
        return s
