import Sun

class Gensun(object):
    """  Generate sun object from astronomical data and location latitude""" 

    def __init__(self):
        pass


    def __call__(self, Rsun,DOY,heureTU,lat):
        s=Sun.Sun()
        s.Rsun=Rsun
        s._set_pos_astro(DOY,heureTU,lat)

        return s
