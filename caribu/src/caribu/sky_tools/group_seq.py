class group_seq(object):
    """  Doc... """ 

    def __init__(self):
        pass


    def __call__(self, tab, indice):
        return (range(0,len(tab[indice][0])-1),)
